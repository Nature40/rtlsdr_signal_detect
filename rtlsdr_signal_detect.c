//
// rtlsdr_signal_detect.c : detect signals within captured rtlsdr sample file
//
// methodology: compute spectral periodogram, look for where power spectral
// density exceeds threshold, count number of transforms exceeding this
// threshold. Observe time, duration, and bandwidth of signal.
//

#include <getopt.h>
#include <limits.h>
#include <liquid/liquid.h>
#include <math.h>
#include <mysql.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define DB_BASE "rteu"
#define DB_RUNS "rteu.runs"
#define DB_TABLE "rteu.signals"

#define BILLION 1000000000L

#define DB_BASE_CREATE "CREATE DATABASE IF NOT EXISTS " DB_BASE ";"

#define DB_RUNS_CREATE                                                         \
    "CREATE TABLE IF NOT EXISTS " DB_RUNS " ("                                 \
    "id int(10) unsigned NOT NULL AUTO_INCREMENT,"                             \
    "device varchar(60) DEFAULT NULL,"                                         \
    "pos_x float DEFAULT NULL,"                                                \
    "pos_y float DEFAULT NULL,"                                                \
    "orientation smallint(5) unsigned DEFAULT NULL,"                           \
    "beam_width smallint(5) unsigned DEFAULT NULL,"                            \
    "gain tinyint(3) unsigned DEFAULT NULL,"                                   \
    "center_freq int(10) unsigned DEFAULT NULL,"                               \
    "freq_range int(10) unsigned DEFAULT NULL,"                                \
    "threshold tinyint(3) unsigned DEFAULT NULL,"                              \
    "fft_bins smallint(5) unsigned DEFAULT NULL,"                              \
    "fft_sample tinyint(3) unsigned DEFAULT NULL,"                             \
    "PRIMARY  KEY(id)"                                                         \
    ") DEFAULT CHARSET = latin1;"

#define DB_TABLE_CREATE                                                        \
    "CREATE TABLE IF NOT EXISTS " DB_TABLE " ("                                \
    "id          bigint(20) unsigned NOT NULL AUTO_INCREMENT,"                 \
    "timestamp   char(30) DEFAULT NULL,"                                       \
    "samples     bigint(20) DEFAULT NULL,"                                     \
    "duration    float DEFAULT NULL,"                                          \
    "signal_freq float DEFAULT NULL,"                                          \
    "signal_bw   float DEFAULT NULL,"                                          \
    "max_signal  float DEFAULT NULL,"                                          \
    "noise       float DEFAULT NULL,"                                          \
    "run         int(10) unsigned NOT NULL,"                                   \
    "PRIMARY KEY (id),"                                                        \
    "KEY run (run),"                                                           \
    "CONSTRAINT signals_ibfk_1 "                                               \
    "FOREIGN KEY (run) REFERENCES " DB_RUNS " (id) ON DELETE CASCADE"          \
    ") DEFAULT CHARSET=latin1;"

float *psd_template;
float *psd;
float *psd_max;
unsigned long long *psd_sample;
int *detect;
int *count;
int *groups;

int timestep;
int nfft;

unsigned long long num_transforms = 0;
int tmp_transforms = 0;

struct timespec t_start;
int keepalive = 300;

int run_id = 0;

int write_to_db = 0;
MYSQL *con;
char *db_host = NULL, *db_user = NULL, *db_pass = NULL;
unsigned int db_port = 0;

struct option longopts[] = {{"sql", no_argument, &write_to_db, 1},
                            {"db_host", required_argument, NULL, 900},
                            {"db_port", required_argument, NULL, 901},
                            {"db_user", required_argument, NULL, 902},
                            {"db_pass", required_argument, NULL, 903},
                            {"db_run_id", required_argument, NULL, 904},
                            {"ll", required_argument, NULL, 905},
                            {"lu", required_argument, NULL, 906},
                            {0, 0, 0, 0}};

// print usage/help message
void usage() {
    fprintf(
        stderr,
        "%s [-b num] [-d name] [-h] [-i file] [-k sec] [-n num] [-r rate] "
        "[-s] "
        "[--sql [--db_host host] --db_user user --db_pass pass] [-t thre]\n",
        __FILE__);
    fprintf(stderr, "  -h                : print this help\n");
    fprintf(stderr, "  -i <file>         : input data filename\n");
    fprintf(stderr, "  -t <threshold>    : detection threshold above psd, "
                    "default: 10 dB\n");
    fprintf(stderr, "  -s                : use STDIN as input\n");
    fprintf(stderr,
            "  -r <rate>         : sampling rate in Hz, default 250000Hz\n");
    fprintf(stderr, "  -b <number>       : number of bins used for fft, "
                    "default is 400\n");
    fprintf(stderr,
            "  -n <number>       : number of samples per fft, default is 50\n");
    fprintf(stderr, "  -k <seconds>      : prints a keep-alive statement every "
                    "<sec> seconds, default is 300\n");
    fprintf(stderr, " --ll <limit>       : set lower limit in seconds for "
                    "signal duration. Shorter signals will not be logged.");
    fprintf(stderr, " --lu <limit>       : set upper limit in seconds for "
                    "signal duration. Longer signals will not be logged.");
    fprintf(stderr, " --sql              : write to database, requires "
                    "--db_user, --db_pass\n");
    fprintf(stderr, " --db_host <host>   : address of SQL server to use, "
                    "default is localhost\n");
    fprintf(
        stderr,
        " --db_port <pass>   : port on which to connect, use 0 if unsure\n");
    fprintf(stderr, " --db_user <user>   : username for SQL server \n");
    fprintf(stderr, " --db_pass <pass>   : matching password\n");
    fprintf(stderr, " --db_run_id <id>   : numeric id of this recording run. "
                    "Used to link it to its metadata in the SQL database");
}

// read samples from file and store into buffer
unsigned int buf_read(FILE *_fid, float complex *_buf, unsigned int _buf_len);

// forward declaration of methods for signal detection
int update_detect(float _threshold);
int update_count();
int update_groups();
int signal_complete(int _group_id);
float get_group_freq(int _group_id);
float get_group_bw(int _group_id);
float get_group_time(int _group_id);
float get_group_max_sig(int _group_id);
float get_group_noise(int _group_id);
unsigned long long get_group_start_time(int _group_id);
int clear_group_count(int _group_id);
int step(float _threshold, unsigned int _sampling_rate, float lowerLimit,
         float upperLimit);
void format_timestamp(const struct timespec _time, char *_buf,
                      const unsigned long _buf_len);
struct timespec time_add(const struct timespec _t1, const struct timespec _t2);
void free_memory();
void open_connection();
void create_tables();

// main program
int main(int argc, char *argv[]) {
    char filename_input[256] = "data/zeidler-2017-08-06/g10_1e_120kHz.dat";
    float threshold = 10.0f;
    char read_from_stdin = 0;
    unsigned long sampling_rate = 250000;
    nfft = 400;
    timestep = 50;
    float lowerLimit = 0, upperLimit = 1;

    // read command-line options
    int dopt;
    while ((dopt = getopt_long(argc, argv, "hi:t:sr:b:n:d:k:", longopts,
                               NULL)) != -1) {
        switch (dopt) {
        case 'h':
            usage();
            return 0;
        case 'i':
            strncpy(filename_input, optarg, 256);
            break;
        case 't':
            threshold = strtof(optarg, NULL);
            break;
        case 's':
            read_from_stdin = 1;
            break;
        case 'r':
            sampling_rate = atoi(optarg);
            break;
        case 'b':
            nfft = atoi(optarg);
            break;
        case 'n':
            timestep = atoi(optarg);
            break;
        case 'k':
            keepalive = atoi(optarg);
            break;
        case 900:
            db_host = optarg;
            break;
        case 901:
            db_port = atoi(optarg);
            break;
        case 902:
            db_user = optarg;
            break;
        case 903:
            db_pass = optarg;
            break;
        case 904:
            run_id = atoi(optarg);
            break;
        case 905:
            lowerLimit = strtof(optarg, NULL);
            break;
        case 906:
            upperLimit = strtof(optarg, NULL);
            break;
        case 0:
            break; // return value of getopt_long() when setting a flag
        default:
            exit(1);
        }
    }

    // reset counters, etc.
    psd_template = (float *)calloc(nfft, sizeof(float));
    psd = (float *)calloc(nfft, sizeof(float));
    psd_max = (float *)calloc(nfft, sizeof(float));
    detect = (int *)calloc(nfft, sizeof(int));
    count = (int *)calloc(nfft, sizeof(int));
    groups = (int *)calloc(nfft, sizeof(int));
    psd_sample = (unsigned long long *)calloc(nfft, sizeof(unsigned long long));
    keepalive *= sampling_rate / timestep;
    keepalive += keepalive % 16;

    // create spectrogram
    spgramcf periodogram =
        spgramcf_create(nfft, LIQUID_WINDOW_HAMMING, nfft / 2, timestep);

    // buffer
    unsigned int buf_len = 64;
    float complex buf[buf_len];

    // DC-blocking filter 1e-3f
    iirfilt_crcf dcblock = iirfilt_crcf_create_dc_blocker(1e-3f);

    // open SQL database
    if (write_to_db != 0) {
        open_connection();
    }

    // open input file
    FILE *fid;
    if (read_from_stdin) {
        fprintf(stderr, "reading from stdin.\n");
        fid = stdin;
    } else {
        fid = fopen(filename_input, "r");
        if (fid == NULL) {
            fprintf(stderr, "error: could not open %s for reading\n",
                    filename_input);
            exit(-1);
        }
    }

    if (write_to_db != 0)
        fprintf(stderr, "Also sending data to SQL Server at %s.\n", db_host);

    fprintf(stderr, "Ignoring signals shorter than %f and longer than %f\n",
            lowerLimit, upperLimit);

    // print row names
    printf("%29s; %7s; %8s; %15s; %11s; %10s; %10s\n", "timestamp", "samples",
           "duration", "signal_freq", "signal_bw", "max_signal", "noise");
    fflush(stdout);

    clock_gettime(CLOCK_REALTIME, &t_start);
    char tbuf[30];
    format_timestamp(t_start, tbuf, 30);
    fprintf(stderr, "Will print timestamp every %i transforms\n", keepalive);
    fprintf(stderr, "%s\n", tbuf);

    // continue processing as long as there are samples in the file
    unsigned long int total_samples = 0;
    num_transforms = 0;
    do {
        // read samples into buffer
        unsigned int r = buf_read(fid, buf, buf_len);
        if (r != buf_len)
            break;

        // apply DC blocking filter
        iirfilt_crcf_execute_block(dcblock, buf, buf_len, buf);

        // accumulate spectrum
        spgramcf_write(periodogram, buf, buf_len);

        // get number of transforms per cycle tmp_transforms
        tmp_transforms = spgramcf_get_num_transforms(periodogram);

        if (tmp_transforms >= 16) {
            // compute power spectral density output
            spgramcf_get_psd(periodogram, psd);

            // compute average template for one second
            if (num_transforms <= sampling_rate / timestep) {
                // set template PSD for relative signal detection
                // Add up all signal strength to derive minimum value
                int i;
                for (i = 0; i < nfft; i++) {
                    if (psd[i] < psd_template[i]) {
                        psd_template[i] = psd[i];
                    }
                }
                memmove(psd_max, psd, nfft * sizeof(float));
            } else {
                // detect differences between current PSD estimate and
                // template
                step(threshold, sampling_rate, lowerLimit, upperLimit);
            }

            // update counters and reset spectrogram object
            num_transforms += spgramcf_get_num_transforms(periodogram);
            spgramcf_reset(periodogram);

            // print keepalive
            if (num_transforms % keepalive == 0) {
                struct timespec now;
                clock_gettime(CLOCK_REALTIME, &now);
                num_transforms = 0;
                t_start = now;
                char tbuf[30];
                char sql_statement[256];
                format_timestamp(now, tbuf, 30);
                printf("%29s; ; ; ; ; \n", tbuf);
                fflush(stdout);
                if (write_to_db != 0) {
                    snprintf(sql_statement, sizeof(sql_statement),
                             "INSERT INTO %s "
                             "(timestamp,samples,duration,signal_freq,signal_"
                             "bw,max_signal,noise,run) "
                             "VALUE(\"%s\",0,0.0,0.0,0.0,0.0,0.0,%i)",
                             DB_TABLE, tbuf, run_id);
                    mysql_query(con, sql_statement);
                    if (*mysql_error(con))
                        fprintf(stderr, "Error while writing to db: %s\n",
                                mysql_error(con));
                }
            }
        }

        // update total sample count
        total_samples += buf_len;

    } while (!feof(fid));

    // close input files
    fclose(fid);
    if (con != NULL)
        mysql_close(con);

    // write accumulated PSD
    spgramcf_destroy(periodogram);
    iirfilt_crcf_destroy(dcblock);

    fprintf(stderr, "total samples in : %lu\n", total_samples);
    fprintf(stderr, "total transforms : %llu\n", num_transforms);

    free_memory();

    return 0;
}

// read samples from file and store into buffer
unsigned int buf_read(FILE *_fid, float complex *_buf, unsigned int _buf_len) {
    int num_read = 0;
    unsigned int i;
    uint8_t buf2[2];

    for (i = 0; i < _buf_len; i++) {
        // try to read 2 samples at a time
        if (fread(buf2, sizeof(uint8_t), 2, _fid) != 2)
            break;
        // convert to float complex type
        float complex x = ((float)(buf2[0]) - 127.0f) +
                          ((float)(buf2[1]) - 127.0f) * _Complex_I;
        // scale resulting samples
        _buf[i] = x * 1e-3f;
        num_read++;
    }
    return num_read;
}

// update detect
int update_detect(float _threshold) {
    int i;
    int total = 0;
    for (i = 0; i < nfft; i++) {
        // relative
        // detect[i] = ((psd[i] - psd_template[i]) > _threshold) ? 1 : 0;
        if ((psd[i] - psd_template[i]) > _threshold) {
            if (psd_sample[i] == 0)
                psd_sample[i] = num_transforms;
            detect[i] = 1; // write matrix for detection
            psd_max[i] = (psd_max[i] > psd[i]) ? psd_max[i]
                                               : psd[i]; // save highest values
        } else {
            detect[i] = 0;
        }
        // absolute
        // detect[i] = (psd[i] > _threshold) ? 1 : 0;
        total += detect[i];
    }
    return total;
}

// update count
int update_count() {
    int i;
    int total = 0;
    for (i = 0; i < nfft; i++) {
        count[i] += detect[i];
        total += count[i];
    }
    return total;
}

// update groups
int update_groups() {
    // replace all non-zero entries with a 1
    int i;
    for (i = 0; i < nfft; i++)
        groups[i] = count[i] > 0;

    // look for adjacent groups and refactor...
    int group_id = 0;
    i = 0;
    while (i < nfft) {
        // find non-zero value
        if (count[i] == 0) {
            i++;
            continue;
        }

        //
        group_id++;
        while (count[i] > 0 && i < nfft) {
            groups[i] = group_id;
            i++;
        }
    }

    // return number of groups
    return group_id;
}

// update signal detection
int signal_complete(int _group_id) {
    int i;
    for (i = 0; i < nfft; i++) {
        if (groups[i] == _group_id && detect[i])
            return 0;
    }
    return 1;
}

// get group center frequency
float get_group_freq(int _group_id) {
    int i, n = 0;
    float fc = 0.0f;
    for (i = 0; i < nfft; i++) {
        if (groups[i] == _group_id) {
            fc += ((float)i / (float)nfft - 0.5f) * count[i];
            n += count[i];
        }
    }
    return fc / (float)n;
}

// get group bandwidth
float get_group_bw(int _group_id) {
    int i, imin = -1, imax = -1;
    for (i = 0; i < nfft; i++) {
        if (groups[i] == _group_id) {
            if (i < imin || imin == -1)
                imin = i;
            if (i > imax || imax == -1)
                imax = i;
        }
    }
    if (imin == -1 || imax == -1)
        return 0.0f;

    return (float)(imax - imin + 1) / (float)nfft;
}

// get group maximum count from group
float get_group_time(int _group_id) {
    int i;
    int max = -1;
    for (i = 0; i < nfft; i++) {
        if (groups[i] == _group_id && count[i] > max)
            max = count[i];
    }
    return max;
}

// get group maximum signal from group
float get_group_max_sig(int _group_id) {
    int i;
    float max = -1000;
    for (i = 0; i < nfft; i++) {
        if (groups[i] == _group_id && psd_max[i] > max) {
            max = psd_max[i];
        }
    }
    return max + 100; // 10*log10(1e10*(max/10) / 1e10*(noise/10))
}

// get noise at max signal
float get_group_noise(int _group_id) {
    int i;
    float noise = 0;
    float max = -1000;
    for (i = 0; i < nfft; i++) {
        if (groups[i] == _group_id && psd_max[i] > max) {
            max = psd_max[i];
            noise = psd_template[i];
        }
    }
    return noise;
}

// get earliest timestamp for given group
unsigned long long get_group_start_time(int _group_id) {
    int i;
    unsigned long long starttime = ULLONG_MAX;
    for (i = 0; i < nfft; i++) {
        if (groups[i] == _group_id && psd_sample[i] < starttime) {
            starttime = psd_sample[i];
        }
    }
    return starttime;
};

// clear count and max for group
int clear_group_count(int _group_id) {
    int i;
    for (i = 0; i < nfft; i++) {
        if (groups[i] == _group_id) {
            count[i] = 0;
            psd_max[i] = -1000;
            psd_sample[i] = 0;
        }
    }
    return 0;
}

// look for signal
int step(float _threshold, unsigned int _sampling_rate, float lowerLimit,
         float upperLimit) {
    update_detect(_threshold);
    update_count();
    int num_groups = update_groups();
    char timestamp[30];
    char sql_statement[256];
    // determine if signal has stopped based on group and detection
    int i;
    float ratio = (float)timestep / (float)_sampling_rate;
    for (i = 1; i <= num_groups; i++) {
        if (signal_complete(i)) {
            // signal started & stopped
            float duration = tmp_transforms * get_group_time(i) * timestep /
                             _sampling_rate; // duration [samples]
            if (duration >= lowerLimit && duration <= upperLimit) {
                float ftime = (float)get_group_start_time(i) * ratio;
                struct timespec tm;
                tm.tv_nsec = (long)(fmodf(ftime, 1) * BILLION);
                tm.tv_sec = (long)ftime;
                tm = time_add(t_start, tm);
                format_timestamp(tm, timestamp, 30);
                float signal_freq =
                    get_group_freq(i) *
                    _sampling_rate; // center frequency estimate (normalized)
                float signal_bw =
                    get_group_bw(i) *
                    _sampling_rate; // bandwidth estimate (normalized)
                float max_signal =
                    get_group_max_sig(i); // maximum signal strength per group
                float noise = get_group_noise(i); // noise level per group
                printf("%29s; %7llu; %8.6f; %15.6f; %11.6f; %10f; %10f\n",
                       timestamp, num_transforms, duration, signal_freq,
                       signal_bw, max_signal, noise);
                fflush(stdout);
                if (write_to_db != 0) {
                    snprintf(sql_statement, sizeof(sql_statement),
                             "INSERT INTO %s "
                             "(timestamp,samples,duration,signal_freq,signal_"
                             "bw, max_signal, noise, run) "
                             "VALUE(\"%s\",%llu,%-10.6f,%9.6f,%9.6f,%f,%f,%i)",
                             DB_TABLE, timestamp, num_transforms, duration,
                             signal_freq, signal_bw, max_signal, noise, run_id);
                    mysql_query(con, sql_statement);
                    if (*mysql_error(con))
                        fprintf(stderr, "Error while writing to db: %s\n",
                                mysql_error(con));
                }
            }
            // reset counters for group
            clear_group_count(i);
        }
    }
    return 0;
}

// pretty-prints _time into _buf
void format_timestamp(const struct timespec _time, char *_buf,
                      const unsigned long _buf_len) {
    char buffer[11];
    const time_t tm = (time_t)_time.tv_sec;
    strftime(_buf, _buf_len, "%F %T", gmtime(&tm));
    sprintf(buffer, ".%09ld", _time.tv_nsec);
    strncat(_buf, buffer, 10);
}

// add 2 timestamps
struct timespec time_add(const struct timespec _t1, const struct timespec _t2) {
    long sec = _t2.tv_sec + _t1.tv_sec;
    long nsec = _t2.tv_nsec + _t1.tv_nsec;
    if (nsec >= BILLION) {
        nsec -= BILLION;
        sec++;
    }
    return (struct timespec){.tv_sec = sec, .tv_nsec = nsec};
}

void free_memory() {
    free(psd_template);
    free(psd);
    free(psd_max);
    free(psd_sample);
    free(detect);
    free(count);
    free(groups);
}

void open_connection() {
    if (db_user == NULL) {
        db_user = "rteu";
        fprintf(stderr, "No database user supplied, trying \"%s\"...\n",
                db_user);
    }
    if (db_pass == NULL) {
        db_pass = "rteu";
        fprintf(stderr, "No database password supplied, trying \"%s\"...\n",
                db_pass);
    }
    if (db_host == NULL) {
        db_host = "localhost";
        fprintf(stderr, "No hostname given, trying \"%s\".\n", db_host);
    }

    con = mysql_init(NULL);
    my_bool reconnect = 1;
    mysql_options(con, MYSQL_OPT_RECONNECT, &reconnect);
    if (con != NULL) {
        if (mysql_real_connect(con, db_host, db_user, db_pass, NULL, db_port,
                               NULL, 0) == NULL) {
            fprintf(stderr, "Error connecting to database: %s\n",
                    mysql_error(con));
            mysql_close(con);
            write_to_db = 0;
        } else {
            create_tables();
        }
    } else {
        fprintf(stderr, "Error connecting to database: %s\n", mysql_error(con));
        write_to_db = 0;
    }
}

void create_tables() {
    char *create_cmds[] = {DB_BASE_CREATE, DB_RUNS_CREATE, DB_TABLE_CREATE};

    for (int i = 0; i < 3; i++) {
        mysql_query(con, create_cmds[i]);
        if (*mysql_error(con)) {
            fprintf(stderr, "Error during create command: %s\n",
                    mysql_error(con));
            write_to_db = 0;
            return;
        }
    }
}