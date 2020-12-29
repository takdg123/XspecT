import os
from subprocess import PIPE, Popen
from threading import Thread
from queue import Queue, Empty  # python 3.x
import time

def enqueue_output(out, queue):
    for line in iter(out.readline, ''):
        queue.put(line)
    out.close()


class XspecInterface:
    
    def __init__(self, logger=print):
        xspec = 'xspec'

        xspec_env = os.environ.copy()
        if 'HEADAS' in xspec_env:
            orig_path = ''
            if 'DYLD_LIBRARY_PATH' in xspec_env:
                orig_path = xspec_env['DYLD_LIBRARY_PATH']
            new_path = os.path.join(xspec_env['HEADAS'], 'lib')
            if new_path not in orig_path:
                if orig_path:
                    xspec_env['DYLD_LIBRARY_PATH'] = new_path + ':' + orig_path
                else:
                    xspec_env['DYLD_LIBRARY_PATH'] = new_path
            # print("NEW DYLD_LIBRARY_PATH = {:}".format(xspec_env['DYLD_LIBRARY_PATH']))

        self.logger = logger
        self.pipe = Popen(xspec, stdin=PIPE, stdout=PIPE, bufsize=1, 
                          universal_newlines=True, env=xspec_env)

        self.q = Queue()
        t = Thread(target=enqueue_output, args=(self.pipe.stdout, self.q))
        t.daemon = True  # thread dies with the program
        t.start()

    def command(self, command_str, add_terminator=True):

        # Add a terminating string the command
        if add_terminator == True:
            if r'\n' in command_str:
                command_str = command_str.replace(r'\n', '; echo Process Complete.;\n')
            else:
                command_str = command_str + "; echo Process Complete.;\n"

        # Issue the command
        self.pipe.stdin.write(command_str)
        # And flush the pipe, just to be safe:
        self.pipe.stdin.flush()

    def response(self):
        # read line without blocking
        collecting = True
        while collecting:
            try:
                line = self.q.get_nowait()
            except Empty:
                collecting = False
            else:
                self.logger(line, end='')

    def response_nowait(self, loggerTK=None):

        # read line without blocking
        collecting = True
        lines = []
        while collecting:
            try:
        
                # Get the next line in the queue that contains the stdout lines
                line = self.q.get_nowait()
            except Empty:
                # collecting = False
                pass
            else:
        
                # Stop collecting when the terminating string is encountered 
                if 'echo Process Complete.' not in line and 'Process Complete.' in line:
                    # self.logger(line, end='')
                    collecting = False
                    break
                else:
        
                    # Don't display the terminating string
                    if 'Process Complete.' in line or r'/bin/echo' in line:
                        pass
                    else:
        
                        # line =  r'<newline> ' + line
                        self.logger(line, end='')
                        lines.append(line)

                        if loggerTK is not None:
                            loggerTK.update()
        return lines

    def capture(self):
        # read line without blocking
        lines = []
        collecting = True
        while collecting:
            try:
                line = self.q.get_nowait()
            except Empty:
                collecting = False
            else:
                lines.append(line)
        return lines

    def capture_nowait(self):

        # read line without blocking
        lines = []
        collecting = True
        while collecting:
            try:
        
                # Get the next line in the queue that contains the stdout lines                
                line = self.q.get_nowait()
            except Empty:
                # print("Waiting")
                pass
            else:
        
                # Stop collecting when the terminating string is encountered 
                if 'echo Process Complete.' not in line and 'Process Complete.' in line:
                    break
                else:
        
                    # Don't capture the terminating string
                    if 'Process Complete.' in line or r'/bin/echo' in line:
                        pass
                    else:
                        lines.append(line.rstrip())


        return lines