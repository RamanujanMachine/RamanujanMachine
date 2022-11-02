'''
'''
from logging import getLogger
from logging.config import fileConfig
import signal
import numpy as np
import os
import sys
from db.config import configuration
from db.lib.pool_handler import WorkerPool

LOGGER_NAME = 'job_logger'
MOD_PATH = 'jobs.job_%s'

def main() -> None:
    os.makedirs(os.path.join(os.getcwd(), 'logs'), exist_ok=True)
    with open('pid.txt', 'w') as pid_file:
        pid_file.writelines([str(os.getpid()), os.linesep])
    worker_pool = WorkerPool(configuration['pool_size'])
    signal.signal(signal.SIGINT, lambda sig, frame: worker_pool.stop())
    results = worker_pool.start({MOD_PATH % name : config for name, config in configuration['jobs_to_run'].items() })
    fileConfig('db/logging.config', defaults={'log_filename': 'main'})

    for module_path, timings in results:
        getLogger(LOGGER_NAME).info('-------------------------------------')
        if timings:
            getLogger(LOGGER_NAME).info(f'module {module_path} running times:')
            getLogger(LOGGER_NAME).info(f'min time: {min(timings)}')
            getLogger(LOGGER_NAME).info(f'max time: {max(timings)}')
            getLogger(LOGGER_NAME).info(f'median time: {np.median(timings)}')
            getLogger(LOGGER_NAME).info(f'average time: {np.average(timings)}')
        else:
            getLogger(LOGGER_NAME).info(f"module {module_path} didn't run! check logs")
        getLogger(LOGGER_NAME).info('-------------------------------------')
        

def stop() -> None:
    print('stopping')
    with open('pid.txt', 'r') as pid_file:
        lines = pid_file.readlines()
    os.kill(int(lines[0].strip()), signal.CTRL_C_EVENT if os.name == 'nt' else signal.SIGINT)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Only commands are start and stop')
        exit(1)

    if sys.argv[1] == 'stop':
        stop()
    elif sys.argv[1] == 'start':
        main()
    else:
        print('Only commands are start and stop')
        exit(1)
