import signal

def signal_handler(signum, frame):
    print 'Signal handler called with signal', signum
    exit(signum)


signals = {
        signal.SIGABRT: 'SIGABRT',
        signal.SIGALRM: 'SIGALRM',
        signal.SIGBUS: 'SIGBUS',
        #signal.SIGCHLD: 'SIGCHLD',
        signal.SIGCONT: 'SIGCONT',
        signal.SIGFPE: 'SIGFPE',
        signal.SIGHUP: 'SIGHUP',
        signal.SIGILL: 'SIGILL',
        signal.SIGINT: 'SIGINT',
        signal.SIGPIPE: 'SIGPIPE',
        signal.SIGPOLL: 'SIGPOLL',
        signal.SIGPROF: 'SIGPROF',
        signal.SIGQUIT: 'SIGQUIT',
        signal.SIGSEGV: 'SIGSEGV',
        signal.SIGSYS: 'SIGSYS',
        signal.SIGTERM: 'SIGTERM',
        signal.SIGTRAP: 'SIGTRAP',
        signal.SIGTSTP: 'SIGTSTP',
        signal.SIGTTIN: 'SIGTTIN',
        signal.SIGTTOU: 'SIGTTOU',
        signal.SIGURG: 'SIGURG',
        signal.SIGUSR1: 'SIGUSR1',
        signal.SIGUSR2: 'SIGUSR2',
        signal.SIGVTALRM: 'SIGVTALRM',
        signal.SIGXCPU: 'SIGXCPU',
        signal.SIGXFSZ: 'SIGXFSZ',
        }

for signum in signals:
    signal.signal(signum, signal_handler)
