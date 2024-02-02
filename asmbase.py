"""
    This is a Python implementation of the ASM2d-N2O model.

    -   Inspired by an implementation of the ASM1 model in the PooPyLab library (https://pypi.org/project/poopylab/).
    -   Convert the ASM2d-N2O model from raw stoichiometric matrix and kinetic equations to a Python class.

    Reference:
        Massara, T.M., Solís, B., Guisasola, A., Katsou, E. and Baeza, J.A., 2018. 
        Development of an ASM2d-N2O model to describe nitrous oxide emissions in municipal WWTPs under dynamic conditions. 
        Chemical Engineering Journal, 335, pp.185-196.
        (https://doi.org/10.1016/j.cej.2017.10.119)
"""


class asm_model(object):


    def __init__(self, ww_temp=20, DO=2):
        '''
        Initialize the ASM2d-N2O model with water temperature and dissolved oxygen.

        Args:
            ww_temp:   wastewater temperature, degC;
            DO:        dissoved oxygen, mg/L

        Return:
            None

        '''

        # wastewater temperature used in the model, degC
        self._temperature = ww_temp
        # mixed liquor bulk dissolved oxygen, mg/L
        self._bulk_DO = DO

        # default kinetic constants AT 20 degree Celcius under ideal conditions
        self._kinetics_20C = {}

        # kinetic parameters AT PROJECT TEMPERATURE
        self._params = {}

        # stoichiometrics
        self._stoichs = {}

        # ASM model components
        self._comps = []

        # temperature difference b/t what's used and baseline (20C), degC
        self._delta_t = self._temperature - 20

        # oxygen transfer coefficient, mg/L-day
        # placeholder for now, let
        #   OUR = 60 mg/L-hr;
        #   Saturation DO = 9 mg/L; and
        #   DO in mixed liquor = 2 mg/L
        self._KLa = 60 * 24 / (9 - 2)

        return None


    def alter_kinetic_20C(self, name, new_val):
        '''
        Alter a model kinetic constant at 20C (baseline temperature).

        Args:
            name:       name of the parameter to be altered (i.e. 'u_max_H');
            new_val:    new parameter numeric value at 20C

        Return:
            None

        '''
        if new_val > 0 and name in self._kinetics_20C.keys():
            self._kinetics_20C[name] = new_val
        else:
            print('ERROR IN ALTERING 20C KINETICS. NO PARAMETER WAS CHANGED')

        return None


    def update(self, ww_temp, DO):
        ''' 
        Update the ASM model with new water temperature and dissolved O2. 

        Args:
            ww_temp:    wastewater temperature, degC;
            DO:         dissolved oxygen, mg/L

        Return:
            None

        '''
        self._temperature = ww_temp
        self._bulk_DO = DO
        self._delta_t = self._temperature - 20.0
        self._set_params()
        self._set_stoichs()
        return None


    def set_KLa(self, kla):
        '''
        Set KLa value.

        Args:
            kla:    new KLa value, mg/m3-day

        Return:
            None
            
        '''
        if kla > 0:
            self._KLa = kla
        else:
            print("ERROR IN USER GIVEN KLa VALUE. KLa NOT CHANGED")

        return None


    def get_params(self):
        '''
        Return the values of the kinetic parameter dictionary.

        '''
        return self._params.copy()


    def get_stoichs(self):
        '''
        Return the values of the stoichiometric dictionary.

        '''
        return self._stoichs.copy()


    def get_all_comps(self):
        '''
        Return a copy of the model components (concentrations).

        '''
        return self._comps[:]


    def get_bulk_DO(self):
        '''
        Return the bulk dissolved O2 concentration, mg/L.

        '''
        return self._bulk_DO


    def _set_ideal_kinetics_20C_to_defaults(self):
        '''
        Set the kinetic params/consts @ 20C to default ideal values.

        '''
        pass


    def _set_params(self):
        '''
        Set the kinetic parameters/constants to the project temperature & DO.

        '''
        pass


    def _set_stoichs(self):
        '''
        Set the stoichiometrics for the model.
        
        '''
        pass


    def _monod(self, term_in_num_denum, term_only_in_denum):
        '''
        Template for Monod kinetics or switches.

        Args:
            term_in_num_denum:      the term in both numerator & denumerator
            term_only_in_denum:     the term only in numerator

        Return:
            float
        '''
        return term_in_num_denum / (term_in_num_denum + term_only_in_denum)


    def _dCdt(self, t, mo_comps, vol, flow, in_comps):
        '''
        Defines dC/dt for the reactor based on mass balance.

        '''
        pass

