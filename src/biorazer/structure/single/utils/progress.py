from progress.bar import Bar
from ipywidgets import IntProgress
from IPython.display import display

class MyBar:
    
    def __init__(self, min, max, description=""):
        self.min = min
        self.max = max
        self.description = description
        try:
            shell = get_ipython().__class__.__name__
            if shell == 'ZMQInteractiveShell' or shell == 'TerminalInteractiveShell':
                # print("Progress bar for Jupyter")
                self.bar = IntProgress(min=self.min, max=self.max, description=self.description)
                display(self.bar)
            else:
                raise EnvironmentError("Not using Jupyter or IPython")
        except NameError:
            # print("Progress bar not using Jupyter")
            self.bar = Bar(self.description, min=self.min, max=self.max)
    
    def next(self, n):
        if isinstance(self.bar, IntProgress):
            self.bar.value += n
        elif isinstance(self.bar, Bar):
            self.bar.next(n)
        else:
            raise ValueError(f"Invalid bar type {type(self.bar)}")
        
    def set_bar_value(self, n):
        if isinstance(self.bar, IntProgress):
            self.bar.value = n
        elif isinstance(self.bar, Bar):
            self.bar.goto(n)
        else:
            raise ValueError(f"Invalid bar type {type(self.bar)}")
        
    def finish(self):
        if isinstance(self.bar, IntProgress):
            self.bar.close()
        elif isinstance(self.bar, Bar):
            self.bar.finish()
        else:
            raise ValueError(f"Invalid bar type {type(self.bar)}")