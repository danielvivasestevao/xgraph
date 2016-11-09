import datetime
from warnings import warn


class ProcessTimer:
    def __init__(self, msg: str=None):
        if msg:
            self.msg = msg
        else:
            self.msg = ""
        self.start_time = None
        self.end_time = None

    def start(self):
        self.start_time = datetime.datetime.now()
        print("[", self.start_time, "]", " START ", self.msg, sep="",
              flush=True)

    def end(self):
        if not self.start_time:
            warn("The timer has not been started yet. Abort")
            return

        self.end_time = datetime.datetime.now()
        print("[", self.end_time, "]", " END   ", self.msg, sep="", flush=True)
        time_delta = self.end_time - self.start_time
        print("[", time_delta, "]", " PROCESSING TIME",
              sep="", flush=True)
        return time_delta
