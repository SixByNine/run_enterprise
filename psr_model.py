import dill
from enterprise.pulsar import Pulsar
from enterprise.signals import signal_base

class PsrModel:
    def __init__(self, model, tim,par,psr=None):
        self.model = model
        self._model = None
        self.tim = tim
        self.par = par
        self.psr = psr
        self.pta = None

    def get_pulsar(self):
        self.psr = Pulsar(self.par, self.tim, drop_t2pulsar=False)

    def get_pta(self):
        if self.model is None:
            self.model = dill.loads(self._model)
        if self.psr is None:
            self.get_pulsar()
        if self.pta is None:
            self.pta = signal_base.PTA(self.model(self.psr))
        return self.pta

    def pack_for_thread(self):
        self._model  = dill.dumps(self.model)
        self.model=None
        self.psr = None
        self.pta = None
