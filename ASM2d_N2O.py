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


import math

from .asmbase import asm_model

class ASM2d_N2O(asm_model):
    '''
    Kinetics and stoichiometrics for ASM2d-N2O model.

    '''

    def __init__(self, ww_temp=20, DO=2):
        '''
        Initialize the ASM2d-N2O model with water temperature and dissolved oxygen.

        Args:
            ww_temp:   wastewater temperature, degC;
            DO:        dissoved oxygen, mg/L

        Return:
            None

        '''
        asm_model.__init__(self)

        self._set_ideal_kinetics_20C_to_defaults()

        # wastewater temperature used in the model, degC
        self._temperature = ww_temp
        # mixed liquor bulk dissolved oxygen, mg/L
        self._bulk_DO = DO

        # temperature gap between project condition and default value, degC
        self._delta_t = self._temperature - 20

        self.update(ww_temp, DO)

        # ASM2d-N2O model components
        self._comps = [0.0] * 24

        # intermediate results of Monod or Inhibition Terms
        self._monods = [1.0] * 52

        # intermediate results of rate expressions
        self._rate_res = [0.0] * 40

        return None
    
    def _set_ideal_kinetics_20C_to_defaults(self):
        '''
        Set ideal kinetics at 20 degree Celcius to default values.
        
        '''
        
        ## Default parameters (conversion factors, stoihiometric coefficients, kinetic parameters) at 20 degree Celcius
        ## definition can be found in the supplementary info of ASM2d-N2O model documentation (https://doi.org/10.1016/j.cej.2017.10.119)
        
        ## Conversion factors
        self._kinetics_20C['i_NSF'] = 0.03
        self._kinetics_20C['i_PSF'] = 0.01
        self._kinetics_20C['i_NSI'] = 0.01
        self._kinetics_20C['i_PSI'] = 0.0
        self._kinetics_20C['i_NXI'] = 0.02
        self._kinetics_20C['i_PXI'] = 0.01
        self._kinetics_20C['i_TSSXI'] = 0.75
        self._kinetics_20C['i_NXS'] = 0.04
        self._kinetics_20C['i_PXS'] = 0.01
        self._kinetics_20C['i_TSSXS'] = 0.75
        self._kinetics_20C['i_NBM'] = 0.07
        self._kinetics_20C['i_PBM'] = 0.02
        self._kinetics_20C['i_TSSBM'] = 0.9

        ## Stoichiometric coefficients
        self._kinetics_20C['Y_H'] = 0.625
        self._kinetics_20C['Y_PHA'] = 0.2
        self._kinetics_20C['Y_PAO'] = 0.625
        self._kinetics_20C['Y_PO4'] = 0.4
        self._kinetics_20C['Y_AOB'] = 0.18
        self._kinetics_20C['Y_NOB'] = 0.08
        self._kinetics_20C['f_SI'] = 0.0
        self._kinetics_20C['f_XI'] = 0.1
        self._kinetics_20C['n_G'] = 1.0

        ## Kinetic parameters
        self._kinetics_20C['K_H'] = 3.0
        self._kinetics_20C['K_O2_H'] = 0.2
        self._kinetics_20C['K_X_H'] = 0.1
        self._kinetics_20C['n_NO3_H'] = 0.6
        self._kinetics_20C['n_NO2_H'] = 0.6
        self._kinetics_20C['K_NO3_H'] = 0.5
        self._kinetics_20C['K_NO2_H'] = 0.5
        self._kinetics_20C['n_fe_H'] = 0.4

        self._kinetics_20C['mu_H'] = 6.0
        self._kinetics_20C['K_O2'] = 0.2
        self._kinetics_20C['K_F'] = 4.0
        self._kinetics_20C['K_NH4'] = 0.05
        self._kinetics_20C['K_P'] = 0.01
        self._kinetics_20C['K_ALK'] = 0.1
        self._kinetics_20C['K_A'] = 4.0
        self._kinetics_20C['K_NO3'] = 0.5
        self._kinetics_20C['K_NO2'] = 0.5
        self._kinetics_20C['n_NO3_D'] = 0.8
        self._kinetics_20C['q_fe'] = 3.0
        self._kinetics_20C['K_fe_H'] = 4.0
        self._kinetics_20C['b_H'] = 0.4
        self._kinetics_20C['mu_H_Den'] = 6.25
        self._kinetics_20C['n_G3'] = 0.16
        self._kinetics_20C['n_G4'] = 0.35
        self._kinetics_20C['n_G5'] = 0.35
        self._kinetics_20C['K_S3'] = 20.0
        self._kinetics_20C['K_S4'] = 20.0
        self._kinetics_20C['K_S5'] = 40.0
        self._kinetics_20C['K_NO2_Den'] = 0.2
        self._kinetics_20C['K_OH4'] = 0.1
        self._kinetics_20C['K_N2O_Den'] = 0.05
        self._kinetics_20C['K_OH3'] = 0.1
        self._kinetics_20C['K_NO_Den'] = 0.05
        self._kinetics_20C['K_OH5'] = 0.1
        self._kinetics_20C['K_I3NO'] = 0.5
        self._kinetics_20C['K_I4NO'] = 0.3
        self._kinetics_20C['K_I5NO'] = 0.075

        self._kinetics_20C['q_PHA'] = 3.0
        self._kinetics_20C['K_A_P'] = 4.0
        self._kinetics_20C['K_ALK_P'] = 0.1
        self._kinetics_20C['q_PP'] = 1.5
        self._kinetics_20C['K_O2_P'] = 0.2
        self._kinetics_20C['K_P_P'] = 0.2
        self._kinetics_20C['K_PHA_P'] = 0.01
        self._kinetics_20C['K_MAX_P'] = 0.34
        self._kinetics_20C['K_PP_P'] = 0.01
        self._kinetics_20C['K_IPP_P'] = 0.02
        self._kinetics_20C['K_PO4_P'] = 0.01
        self._kinetics_20C['n_NO3_P'] = 0.6
        self._kinetics_20C['n_NO2_P'] = 0.6
        self._kinetics_20C['K_NO3_P'] = 0.5
        self._kinetics_20C['K_NO2_P'] = 0.5
        self._kinetics_20C['mu_PAO'] = 1.0
        self._kinetics_20C['b_PAO'] = 0.2
        self._kinetics_20C['b_PP'] = 0.2
        self._kinetics_20C['b_PHA'] = 0.2

        self._kinetics_20C['mu_AOB_HAO'] = 0.78
        self._kinetics_20C['q_AOB_AMO'] = 5.2008
        self._kinetics_20C['K_O2_AOB1'] = 1.0
        self._kinetics_20C['K_NH4_AOB'] = 0.2
        self._kinetics_20C['K_O2_AOB2'] = 0.6
        self._kinetics_20C['K_NH2OH_AOB'] = 0.9
        self._kinetics_20C['q_AOB_HAO'] = 5.2008
        self._kinetics_20C['K_NO_AOB_HAO'] = 0.0003
        self._kinetics_20C['q_AOB_N2O_NN'] = 0.0078
        self._kinetics_20C['K_NO_AOB_NN'] = 0.008
        self._kinetics_20C['K_O2_AOB_ND'] = 0.5
        self._kinetics_20C['K_I_O2_AOB'] = 0.8
        self._kinetics_20C['K_HNO2_AOB'] = 0.004
        self._kinetics_20C['q_AOB_N2O_ND'] = 1.3008
        self._kinetics_20C['K_ALK_AOB'] = 0.1
        self._kinetics_20C['K_P_AOB'] = 0.01
        self._kinetics_20C['mu_NOB'] = 0.78 # 1.02
        self._kinetics_20C['K_O2_NOB'] = 1.2 # 1.75
        self._kinetics_20C['K_ALK_NOB'] = 0.1
        self._kinetics_20C['K_NO2_NOB'] = 0.5
        self._kinetics_20C['K_P_NOB'] = 0.01
        self._kinetics_20C['b_AOB'] = 0.096
        self._kinetics_20C['b_NOB'] = 0.096 # 0.17

        self._kinetics_20C['k_PRE'] = 1.0
        self._kinetics_20C['k_RED'] = 0.6
        self._kinetics_20C['K_ALK_PR'] = 0.5

        return None
    
    def _set_params(self):
        '''
        Set the kinetic parameters/constants at project temperature.
        TODO: update the parameters based on the temperature by Arrhenius equation.

        '''

        ## Conversion factors
        self._params['i_NSF'] = self._kinetics_20C['i_NSF']
        self._params['i_PSF'] = self._kinetics_20C['i_PSF']
        self._params['i_NSI'] = self._kinetics_20C['i_NSI']
        self._params['i_PSI'] = self._kinetics_20C['i_PSI']
        self._params['i_NXI'] = self._kinetics_20C['i_NXI']
        self._params['i_PXI'] = self._kinetics_20C['i_PXI']
        self._params['i_TSSXI'] = self._kinetics_20C['i_TSSXI']
        self._params['i_NXS'] = self._kinetics_20C['i_NXS']
        self._params['i_PXS'] = self._kinetics_20C['i_PXS']
        self._params['i_TSSXS'] = self._kinetics_20C['i_TSSXS']
        self._params['i_NBM'] = self._kinetics_20C['i_NBM']
        self._params['i_PBM'] = self._kinetics_20C['i_PBM']
        self._params['i_TSSBM'] = self._kinetics_20C['i_TSSBM']

        ## Stoichiometric coefficients
        self._params['Y_H'] = self._kinetics_20C['Y_H']
        self._params['Y_PHA'] = self._kinetics_20C['Y_PHA']
        self._params['Y_PAO'] = self._kinetics_20C['Y_PAO']
        self._params['Y_PO4'] = self._kinetics_20C['Y_PO4']
        self._params['Y_AOB'] = self._kinetics_20C['Y_AOB']
        self._params['Y_NOB'] = self._kinetics_20C['Y_NOB']
        self._params['f_SI'] = self._kinetics_20C['f_SI']
        self._params['f_XI'] = self._kinetics_20C['f_XI']
        self._params['n_G'] = self._kinetics_20C['n_G']

        ## Kinetic parameters
        self._params['K_H'] = self._kinetics_20C['K_H']
        self._params['K_O2_H'] = self._kinetics_20C['K_O2_H']
        self._params['K_X_H'] = self._kinetics_20C['K_X_H']
        self._params['n_NO3_H'] = self._kinetics_20C['n_NO3_H']
        self._params['n_NO2_H'] = self._kinetics_20C['n_NO2_H']
        self._params['K_NO3_H'] = self._kinetics_20C['K_NO3_H']
        self._params['K_NO2_H'] = self._kinetics_20C['K_NO2_H']
        self._params['n_fe_H'] = self._kinetics_20C['n_fe_H']

        self._params['mu_H'] = self._kinetics_20C['mu_H']
        self._params['K_O2'] = self._kinetics_20C['K_O2']
        self._params['K_F'] = self._kinetics_20C['K_F']
        self._params['K_NH4'] = self._kinetics_20C['K_NH4']
        self._params['K_P'] = self._kinetics_20C['K_P']
        self._params['K_ALK'] = self._kinetics_20C['K_ALK']
        self._params['K_A'] = self._kinetics_20C['K_A']
        self._params['K_NO3'] = self._kinetics_20C['K_NO3']
        self._params['K_NO2'] = self._kinetics_20C['K_NO2']
        self._params['n_NO3_D'] = self._kinetics_20C['n_NO3_D']
        self._params['q_fe'] = self._kinetics_20C['q_fe']
        self._params['K_fe_H'] = self._kinetics_20C['K_fe_H']
        self._params['b_H'] = self._kinetics_20C['b_H']
        self._params['mu_H_Den'] = self._kinetics_20C['mu_H_Den']
        self._params['n_G3'] = self._kinetics_20C['n_G3']
        self._params['n_G4'] = self._kinetics_20C['n_G4']
        self._params['n_G5'] = self._kinetics_20C['n_G5']
        self._params['K_S3'] = self._kinetics_20C['K_S3']
        self._params['K_S4'] = self._kinetics_20C['K_S4']
        self._params['K_S5'] = self._kinetics_20C['K_S5']
        self._params['K_NO2_Den'] = self._kinetics_20C['K_NO2_Den']
        self._params['K_OH4'] = self._kinetics_20C['K_OH4']
        self._params['K_N2O_Den'] = self._kinetics_20C['K_N2O_Den']
        self._params['K_OH3'] = self._kinetics_20C['K_OH3']
        self._params['K_NO_Den'] = self._kinetics_20C['K_NO_Den']
        self._params['K_OH5'] = self._kinetics_20C['K_OH5']
        self._params['K_I3NO'] = self._kinetics_20C['K_I3NO']
        self._params['K_I4NO'] = self._kinetics_20C['K_I4NO']
        self._params['K_I5NO'] = self._kinetics_20C['K_I5NO']

        self._params['q_PHA'] = self._kinetics_20C['q_PHA']
        self._params['K_A_P'] = self._kinetics_20C['K_A_P']
        self._params['K_ALK_P'] = self._kinetics_20C['K_ALK_P']
        self._params['q_PP'] = self._kinetics_20C['q_PP']
        self._params['K_O2_P'] = self._kinetics_20C['K_O2_P']
        self._params['K_P_P'] = self._kinetics_20C['K_P_P']
        self._params['K_PHA_P'] = self._kinetics_20C['K_PHA_P']
        self._params['K_MAX_P'] = self._kinetics_20C['K_MAX_P']
        self._params['K_PP_P'] = self._kinetics_20C['K_PP_P']
        self._params['K_IPP_P'] = self._kinetics_20C['K_IPP_P']
        self._params['K_PO4_P'] = self._kinetics_20C['K_PO4_P']
        self._params['n_NO3_P'] = self._kinetics_20C['n_NO3_P']
        self._params['n_NO2_P'] = self._kinetics_20C['n_NO2_P']
        self._params['K_NO3_P'] = self._kinetics_20C['K_NO3_P']
        self._params['K_NO2_P'] = self._kinetics_20C['K_NO2_P']
        self._params['mu_PAO'] = self._kinetics_20C['mu_PAO']
        self._params['b_PAO'] = self._kinetics_20C['b_PAO']
        self._params['b_PP'] = self._kinetics_20C['b_PP']
        self._params['b_PHA'] = self._kinetics_20C['b_PHA']

        self._params['mu_AOB_HAO'] = self._kinetics_20C['mu_AOB_HAO']
        self._params['q_AOB_AMO'] = self._kinetics_20C['q_AOB_AMO']
        self._params['K_O2_AOB1'] = self._kinetics_20C['K_O2_AOB1']
        self._params['K_NH4_AOB'] = self._kinetics_20C['K_NH4_AOB']
        self._params['K_O2_AOB2'] = self._kinetics_20C['K_O2_AOB2']
        self._params['K_NH2OH_AOB'] = self._kinetics_20C['K_NH2OH_AOB']
        self._params['q_AOB_HAO'] = self._kinetics_20C['q_AOB_HAO']
        self._params['K_NO_AOB_HAO'] = self._kinetics_20C['K_NO_AOB_HAO']
        self._params['q_AOB_N2O_NN'] = self._kinetics_20C['q_AOB_N2O_NN']
        self._params['K_NO_AOB_NN'] = self._kinetics_20C['K_NO_AOB_NN']
        self._params['K_O2_AOB_ND'] = self._kinetics_20C['K_O2_AOB_ND']
        self._params['K_I_O2_AOB'] = self._kinetics_20C['K_I_O2_AOB']
        self._params['K_HNO2_AOB'] = self._kinetics_20C['K_HNO2_AOB']
        self._params['q_AOB_N2O_ND'] = self._kinetics_20C['q_AOB_N2O_ND']
        self._params['K_ALK_AOB'] = self._kinetics_20C['K_ALK_AOB']
        self._params['K_P_AOB'] = self._kinetics_20C['K_P_AOB']
        self._params['mu_NOB'] = self._kinetics_20C['mu_NOB']
        self._params['K_O2_NOB'] = self._kinetics_20C['K_O2_NOB']
        self._params['K_ALK_NOB'] = self._kinetics_20C['K_ALK_NOB']
        self._params['K_NO2_NOB'] = self._kinetics_20C['K_NO2_NOB']
        self._params['K_P_NOB'] = self._kinetics_20C['K_P_NOB']
        self._params['b_AOB'] = self._kinetics_20C['b_AOB']
        self._params['b_NOB'] = self._kinetics_20C['b_NOB']

        self._params['k_PRE'] = self._kinetics_20C['k_PRE']
        self._params['k_RED'] = self._kinetics_20C['k_RED']
        self._params['K_ALK_PR'] = self._kinetics_20C['K_ALK_PR']
        
        return None

    def _set_stoichs(self):
        '''
        Set the stoichiometrics for the model.

        '''

        ## Stoichiometric coefficients
        ## definition can be found in the Stoichiometric matrix of ASM2d-N2O model documentation (https://doi.org/10.1016/j.cej.2017.10.119)
        ## _stoichs['x_y'] ==> x is process index, and y is component index
        
        self._stoichs['1_2'] = 1.0 - self._params['f_SI']
        self._stoichs['1_4'] = self._params['i_NXS'] - (1.0 - self._params['f_SI']) * self._params['i_NSF']
        self._stoichs['1_10'] = self._params['i_PXS'] - (1.0 - self._params['f_SI']) * self._params['i_PSF']
        self._stoichs['1_11'] = self._params['f_SI']
        self._stoichs['1_15'] = -1.0
        self._stoichs['1_12'] = (1 / 14) * self._stoichs['1_4'] + (-1.5 / 31) * self._stoichs['1_10']
        self._stoichs['1_22'] = self._params['i_TSSXS'] * self._stoichs['1_15']

        self._stoichs['2_2'] = 1.0 - self._params['f_SI']
        self._stoichs['2_4'] = self._params['i_NXS'] - (1.0 - self._params['f_SI']) * self._params['i_NSF']
        self._stoichs['2_10'] = self._params['i_PXS'] - (1.0 - self._params['f_SI']) * self._params['i_PSF']
        self._stoichs['2_11'] = self._params['f_SI']
        self._stoichs['2_15'] = -1.0
        self._stoichs['2_12'] = (1 / 14) * self._stoichs['2_4'] + (-1.5 / 31) * self._stoichs['2_10']
        self._stoichs['2_22'] = self._params['i_TSSXS'] * self._stoichs['2_15']

        self._stoichs['3_2'] = 1.0 - self._params['f_SI']
        self._stoichs['3_4'] = self._params['i_NXS'] - (1.0 - self._params['f_SI']) * self._params['i_NSF']
        self._stoichs['3_10'] = self._params['i_PXS'] - (1.0 - self._params['f_SI']) * self._params['i_PSF']
        self._stoichs['3_11'] = self._params['f_SI']
        self._stoichs['3_15'] = -1.0
        self._stoichs['3_12'] = (1 / 14) * self._stoichs['3_4'] + (-1.5 / 31) * self._stoichs['3_10']
        self._stoichs['3_22'] = self._params['i_TSSXS'] * self._stoichs['3_15']

        self._stoichs['4_2'] = 1.0 - self._params['f_SI']
        self._stoichs['4_4'] = self._params['i_NXS'] - (1.0 - self._params['f_SI']) * self._params['i_NSF']
        self._stoichs['4_10'] = self._params['i_PXS'] - (1.0 - self._params['f_SI']) * self._params['i_PSF']
        self._stoichs['4_11'] = self._params['f_SI']
        self._stoichs['4_15'] = -1.0
        self._stoichs['4_12'] = (1 / 14) * self._stoichs['4_4'] + (-1.5 / 31) * self._stoichs['4_10']
        self._stoichs['4_22'] = self._params['i_TSSXS'] * self._stoichs['4_15']

        self._stoichs['5_1'] = 1.0 - (1.0 / self._params['Y_H'])
        self._stoichs['5_2'] = -1.0 / self._params['Y_H']
        self._stoichs['5_4'] = self._params['i_NSF'] / self._params['Y_H'] - self._params['i_NBM']
        self._stoichs['5_10'] = self._params['i_PSF'] / self._params['Y_H'] - self._params['i_PBM']
        self._stoichs['5_16'] = 1.0
        self._stoichs['5_12'] = (1 / 14) * self._stoichs['5_4'] + (-1.5 / 31) * self._stoichs['5_10']
        self._stoichs['5_22'] = self._params['i_TSSBM'] * self._stoichs['5_16']

        self._stoichs['6_1'] = 1.0 - (1.0 / self._params['Y_H'])
        self._stoichs['6_3'] = -1.0 / self._params['Y_H']
        self._stoichs['6_4'] = 0.0 - self._params['i_NBM']
        self._stoichs['6_10'] = 0.0 - self._params['i_PBM']
        self._stoichs['6_16'] = 1.0
        self._stoichs['6_12'] = (-1 / 64) * self._stoichs['6_3'] + (1 / 14) * self._stoichs['6_4'] + (-1.5 / 31) * self._stoichs['6_10']
        self._stoichs['6_22'] = self._params['i_TSSBM'] * self._stoichs['6_16']

        self._stoichs['7_2'] = -1.0 / (self._params['Y_H'] * self._params['n_G'])
        self._stoichs['7_4'] = self._params['i_NSF'] / self._params['Y_H'] - self._params['i_NBM']
        self._stoichs['7_8'] = (1.0 - self._params['Y_H'] * self._params['n_G']) / ((8 / 7) * self._params['Y_H'] * self._params['n_G'])
        self._stoichs['7_9'] = 0.0 - (1.0 - self._params['Y_H'] * self._params['n_G']) / ((8 / 7) * self._params['Y_H'] * self._params['n_G'])
        self._stoichs['7_10'] = self._params['i_PSF'] / self._params['Y_H'] - self._params['i_PBM']
        self._stoichs['7_16'] = 1.0
        self._stoichs['7_12'] = (1 / 14) * self._stoichs['7_4'] + (-1 / 14) * self._stoichs['7_8'] + (-1 / 14) * self._stoichs['7_9'] + (-1.5 / 31) * self._stoichs['7_10']
        self._stoichs['7_22'] = self._params['i_TSSBM'] * self._stoichs['7_16']

        self._stoichs['8_2'] = -1.0 / (self._params['Y_H'] * self._params['n_G'])
        self._stoichs['8_4'] = self._params['i_NSF'] / self._params['Y_H'] - self._params['i_NBM']
        self._stoichs['8_7'] = (1.0 - self._params['Y_H'] * self._params['n_G']) / ((4 / 7) * self._params['Y_H'] * self._params['n_G'])
        self._stoichs['8_8'] = 0.0 - (1.0 - self._params['Y_H'] * self._params['n_G']) / ((4 / 7) * self._params['Y_H'] * self._params['n_G'])
        self._stoichs['8_10'] = self._params['i_PSF'] / self._params['Y_H'] - self._params['i_PBM']
        self._stoichs['8_16'] = 1.0
        self._stoichs['8_12'] = (1 / 14) * self._stoichs['8_4'] + (-1 / 14) * self._stoichs['8_8'] + (-1.5 / 31) * self._stoichs['8_10']
        self._stoichs['8_22'] = self._params['i_TSSBM'] * self._stoichs['8_16']

        self._stoichs['9_2'] = -1.0 / (self._params['Y_H'] * self._params['n_G'])
        self._stoichs['9_4'] = self._params['i_NSF'] / self._params['Y_H'] - self._params['i_NBM']
        self._stoichs['9_6'] = (1.0 - self._params['Y_H'] * self._params['n_G']) / ((4 / 7) * self._params['Y_H'] * self._params['n_G'])
        self._stoichs['9_7'] = 0.0 - (1.0 - self._params['Y_H'] * self._params['n_G']) / ((4 / 7) * self._params['Y_H'] * self._params['n_G'])
        self._stoichs['9_10'] = self._params['i_PSF'] / self._params['Y_H'] - self._params['i_PBM']
        self._stoichs['9_16'] = 1.0
        self._stoichs['9_12'] = (1 / 14) * self._stoichs['9_4']+ (-1.5 / 31) * self._stoichs['9_10']
        self._stoichs['9_22'] = self._params['i_TSSBM'] * self._stoichs['9_16']

        self._stoichs['10_2'] = -1.0 / (self._params['Y_H'] * self._params['n_G'])
        self._stoichs['10_4'] = self._params['i_NSF'] / self._params['Y_H'] - self._params['i_NBM']
        self._stoichs['10_6'] = 0.0 - (1.0 - self._params['Y_H'] * self._params['n_G']) / ((4 / 7) * self._params['Y_H'] * self._params['n_G'])
        self._stoichs['10_10'] = self._params['i_PSF'] / self._params['Y_H'] - self._params['i_PBM']
        self._stoichs['10_13'] = (1.0 - self._params['Y_H'] * self._params['n_G']) / ((4 / 7) * self._params['Y_H'] * self._params['n_G'])
        self._stoichs['10_16'] = 1.0
        self._stoichs['10_12'] = (1 / 14) * self._stoichs['10_4'] + (-1.5 / 31) * self._stoichs['10_10']
        self._stoichs['10_22'] = self._params['i_TSSBM'] * self._stoichs['10_16']

        self._stoichs['11_3'] = -1.0 / (self._params['Y_H'] * self._params['n_G'])
        self._stoichs['11_4'] = 0.0 - self._params['i_NBM']
        self._stoichs['11_8'] = (1.0 - self._params['Y_H'] * self._params['n_G']) / ((8 / 7) * self._params['Y_H'] * self._params['n_G'])
        self._stoichs['11_9'] = 0.0 - (1.0 - self._params['Y_H'] * self._params['n_G']) / ((8 / 7) * self._params['Y_H'] * self._params['n_G'])
        self._stoichs['11_10'] = 0.0 - self._params['i_PBM']
        self._stoichs['11_16'] = 1.0
        self._stoichs['11_12'] = (-1 / 64) * self._stoichs['11_3'] + (1 / 14) * self._stoichs['11_4'] + (-1 / 14) * self._stoichs['11_8'] + (-1 / 14) * self._stoichs['11_9'] + (-1.5 / 31) * self._stoichs['11_10']
        self._stoichs['11_22'] = self._params['i_TSSBM'] * self._stoichs['11_16']

        self._stoichs['12_3'] = -1.0 / (self._params['Y_H'] * self._params['n_G'])
        self._stoichs['12_4'] = 0.0 - self._params['i_NBM']
        self._stoichs['12_7'] = (1.0 - self._params['Y_H'] * self._params['n_G']) / ((4 / 7) * self._params['Y_H'] * self._params['n_G'])
        self._stoichs['12_8'] = 0.0 - (1.0 - self._params['Y_H'] * self._params['n_G']) / ((4 / 7) * self._params['Y_H'] * self._params['n_G'])
        self._stoichs['12_10'] = 0.0 - self._params['i_PBM']
        self._stoichs['12_16'] = 1.0
        self._stoichs['12_12'] = (-1 / 64) * self._stoichs['12_3'] + (1 / 14) * self._stoichs['12_4'] + (-1 / 14) * self._stoichs['12_8'] + (-1.5 / 31) * self._stoichs['12_10']
        self._stoichs['12_22'] = self._params['i_TSSBM'] * self._stoichs['12_16']

        self._stoichs['13_3'] = -1.0 / (self._params['Y_H'] * self._params['n_G'])
        self._stoichs['13_4'] = 0.0 - self._params['i_NBM']
        self._stoichs['13_6'] = (1.0 - self._params['Y_H'] * self._params['n_G']) / ((4 / 7) * self._params['Y_H'] * self._params['n_G'])
        self._stoichs['13_7'] = 0.0 - (1.0 - self._params['Y_H'] * self._params['n_G']) / ((4 / 7) * self._params['Y_H'] * self._params['n_G'])
        self._stoichs['13_10'] = 0.0 - self._params['i_PBM']
        self._stoichs['13_16'] = 1.0
        self._stoichs['13_12'] = (-1 / 64) * self._stoichs['13_3'] + (1 / 14) * self._stoichs['13_4'] + (-1.5 / 31) * self._stoichs['13_10']
        self._stoichs['13_22'] = self._params['i_TSSBM'] * self._stoichs['13_16']

        self._stoichs['14_3'] = -1.0 / (self._params['Y_H'] * self._params['n_G'])
        self._stoichs['14_4'] = 0.0 - self._params['i_NBM']
        self._stoichs['14_6'] = 0.0 - (1.0 - self._params['Y_H'] * self._params['n_G']) / ((4 / 7) * self._params['Y_H'] * self._params['n_G'])
        self._stoichs['14_10'] = 0.0 - self._params['i_PBM']
        self._stoichs['14_13'] = (1.0 - self._params['Y_H'] * self._params['n_G']) / ((4 / 7) * self._params['Y_H'] * self._params['n_G'])
        self._stoichs['14_16'] = 1.0
        self._stoichs['14_12'] = (-1 / 64) * self._stoichs['14_3'] + (1 / 14) * self._stoichs['14_4'] + (-1.5 / 31) * self._stoichs['14_10']
        self._stoichs['14_22'] = self._params['i_TSSBM'] * self._stoichs['14_16']

        self._stoichs['15_2'] = -1.0
        self._stoichs['15_3'] = 1.0
        self._stoichs['15_4'] = self._params['i_NSF']
        self._stoichs['15_10'] = self._params['i_PSF']
        self._stoichs['15_12'] = (-1 / 64) * self._stoichs['15_3'] + (1 / 14) * self._stoichs['15_4'] + (-1.5 / 31) * self._stoichs['15_10']

        self._stoichs['16_4'] = self._params['i_NBM'] - self._params['i_NXI'] * self._params['f_XI'] - (1.0 - self._params['f_XI']) * self._params['i_NXS']
        self._stoichs['16_10'] = self._params['i_PBM'] - self._params['i_PXI'] * self._params['f_XI'] - (1.0 - self._params['f_XI']) * self._params['i_PXS']
        self._stoichs['16_14'] = self._params['f_XI']
        self._stoichs['16_15'] = 1.0 - self._params['f_XI']
        self._stoichs['16_16'] = -1.0
        self._stoichs['16_12'] = (1 / 14) * self._stoichs['16_4'] + (-1.5 / 31) * self._stoichs['16_10']
        self._stoichs['16_22'] = self._params['i_TSSXI'] * self._stoichs['16_14'] + self._params['i_TSSXS'] * self._stoichs['16_15'] + self._params['i_TSSBM'] * self._stoichs['16_16']

        self._stoichs['17_3'] = -1.0
        self._stoichs['17_10'] = self._params['Y_PO4']
        self._stoichs['17_18'] = 0.0 - self._params['Y_PO4']
        self._stoichs['17_19'] = 1.0
        self._stoichs['17_12'] = (-1 / 64) * self._stoichs['17_3'] + (-1.5 / 31) * self._stoichs['17_10'] + (-1 / 31) * self._stoichs['17_18']
        self._stoichs['17_22'] = 3.23 * self._stoichs['17_18'] + 0.6 * self._stoichs['17_19']

        self._stoichs['18_1'] = 0.0 - self._params['Y_PHA']
        self._stoichs['18_10'] = -1.0
        self._stoichs['18_18'] = 1.0
        self._stoichs['18_19'] = 0.0 - self._params['Y_PHA']
        self._stoichs['18_12'] = (-1.5 / 31) * self._stoichs['18_10'] + (-1 / 31) * self._stoichs['18_18']
        self._stoichs['18_22'] = 3.23 * self._stoichs['18_18'] + 0.6 * self._stoichs['18_19']

        self._stoichs['19_8'] = self._params['Y_PHA'] / (8 / 7)
        self._stoichs['19_9'] = 0.0 - self._params['Y_PHA'] / (8 / 7)
        self._stoichs['19_10'] = -1.0
        self._stoichs['19_18'] = 1.0
        self._stoichs['19_19'] = 0.0 - self._params['Y_PHA']
        self._stoichs['19_12'] = (-1 / 14) * self._stoichs['19_8'] + (-1 / 14) * self._stoichs['19_9'] + (-1.5 / 31) * self._stoichs['19_10'] + (-1 / 31) * self._stoichs['19_18']
        self._stoichs['19_22'] = 3.23 * self._stoichs['19_18'] + 0.6 * self._stoichs['19_19']

        self._stoichs['20_7'] = self._params['Y_PHA'] / (4 / 7)
        self._stoichs['20_8'] = 0.0 - self._params['Y_PHA'] / (4 / 7)
        self._stoichs['20_10'] = -1.0
        self._stoichs['20_18'] = 1.0
        self._stoichs['20_19'] = 0.0 - self._params['Y_PHA']
        self._stoichs['20_12'] = (-1 / 14) * self._stoichs['20_8'] + (-1.5 / 31) * self._stoichs['20_10'] + (-1 / 31) * self._stoichs['20_18']
        self._stoichs['20_22'] = 3.23 * self._stoichs['20_18'] + 0.6 * self._stoichs['20_19']

        self._stoichs['21_6'] = self._params['Y_PHA'] / (4 / 7)
        self._stoichs['21_7'] = 0.0 - self._params['Y_PHA'] / (4 / 7)
        self._stoichs['21_10'] = -1.0
        self._stoichs['21_18'] = 1.0
        self._stoichs['21_19'] = 0.0 - self._params['Y_PHA']
        self._stoichs['21_12'] = (-1.5 / 31) * self._stoichs['21_10'] + (-1 / 31) * self._stoichs['21_18']
        self._stoichs['21_22'] = 3.23 * self._stoichs['21_18'] + 0.6 * self._stoichs['21_19']

        self._stoichs['22_6'] = 0.0 - self._params['Y_PHA'] / (4 / 7)
        self._stoichs['22_10'] = -1.0
        self._stoichs['22_13'] = self._params['Y_PHA'] / (4 / 7)
        self._stoichs['22_18'] = 1.0
        self._stoichs['22_19'] = 0.0 - self._params['Y_PHA']
        self._stoichs['22_12'] = (-1.5 / 31) * self._stoichs['22_10'] + (-1 / 31) * self._stoichs['22_18']
        self._stoichs['22_22'] = 3.23 * self._stoichs['22_18'] + 0.6 * self._stoichs['22_19']

        self._stoichs['23_1'] = 1.0 - 1.0 / self._params['Y_PAO']
        self._stoichs['23_4'] = 0.0 - self._params['i_NBM']
        self._stoichs['23_10'] = 0.0 - self._params['i_PBM']
        self._stoichs['23_17'] = 1.0
        self._stoichs['23_19'] = -1.0 / self._params['Y_PAO']
        self._stoichs['23_12'] = (1 / 14) * self._stoichs['23_4'] + (-1.5 / 31) * self._stoichs['23_10']
        self._stoichs['23_22'] = self._params['i_TSSBM'] * self._stoichs['23_17'] + 0.6 * self._stoichs['23_19']

        self._stoichs['24_4'] = 0.0 - self._params['i_NBM']
        self._stoichs['24_8'] = (1.0 - self._params['Y_PAO'] * self._params['n_G']) / ((8 / 7) * self._params['Y_PAO'] * self._params['n_G'])
        self._stoichs['24_9'] = 0.0 - (1.0 - self._params['Y_PAO'] * self._params['n_G']) / ((8 / 7) * self._params['Y_PAO'] * self._params['n_G'])
        self._stoichs['24_10'] = 0.0 - self._params['i_PBM']
        self._stoichs['24_17'] = 1.0
        self._stoichs['24_19'] = -1.0 / self._params['Y_PAO']
        self._stoichs['24_12'] = (1 / 14) * self._stoichs['24_4'] + (-1 / 14) * self._stoichs['24_8'] + (-1 / 14) * self._stoichs['24_9'] + (-1.5 / 31) * self._stoichs['24_10']
        self._stoichs['24_22'] = self._params['i_TSSBM'] * self._stoichs['24_17'] + 0.6 * self._stoichs['24_19']

        self._stoichs['25_4'] = 0.0 - self._params['i_NBM']
        self._stoichs['25_7'] = (1.0 - self._params['Y_PAO'] * self._params['n_G']) / ((4 / 7) * self._params['Y_PAO'] * self._params['n_G'])
        self._stoichs['25_8'] = 0.0 - (1.0 - self._params['Y_PAO'] * self._params['n_G']) / ((4 / 7) * self._params['Y_PAO'] * self._params['n_G'])
        self._stoichs['25_10'] = 0.0 - self._params['i_PBM']
        self._stoichs['25_17'] = 1.0
        self._stoichs['25_19'] = -1.0 / self._params['Y_PAO']
        self._stoichs['25_12'] = (1 / 14) * self._stoichs['25_4'] + (-1 / 14) * self._stoichs['25_8'] + (-1.5 / 31) * self._stoichs['25_10']
        self._stoichs['25_22'] = self._params['i_TSSBM'] * self._stoichs['25_17'] + 0.6 * self._stoichs['25_19']

        self._stoichs['26_4'] = 0.0 - self._params['i_NBM']
        self._stoichs['26_6'] = (1.0 - self._params['Y_PAO'] * self._params['n_G']) / ((4 / 7) * self._params['Y_PAO'] * self._params['n_G'])
        self._stoichs['26_7'] = 0.0 - (1.0 - self._params['Y_PAO'] * self._params['n_G']) / ((4 / 7) * self._params['Y_PAO'] * self._params['n_G'])
        self._stoichs['26_10'] = 0.0 - self._params['i_PBM']
        self._stoichs['26_17'] = 1.0
        self._stoichs['26_19'] = -1.0 / self._params['Y_PAO']
        self._stoichs['26_12'] = (1 / 14) * self._stoichs['26_4'] + (-1.5 / 31) * self._stoichs['26_10']
        self._stoichs['26_22'] = self._params['i_TSSBM'] * self._stoichs['26_17'] + 0.6 * self._stoichs['26_19']

        self._stoichs['27_4'] = 0.0 - self._params['i_NBM']
        self._stoichs['27_6'] = 0.0 - (1.0 - self._params['Y_PAO'] * self._params['n_G']) / ((4 / 7) * self._params['Y_PAO'] * self._params['n_G'])
        self._stoichs['27_10'] = 0.0 - self._params['i_PBM']
        self._stoichs['27_13'] = (1.0 - self._params['Y_PAO'] * self._params['n_G']) / ((4 / 7) * self._params['Y_PAO'] * self._params['n_G'])
        self._stoichs['27_17'] = 1.0
        self._stoichs['27_19'] = -1.0 / self._params['Y_PAO']
        self._stoichs['27_12'] = (1 / 14) * self._stoichs['27_4'] + (-1.5 / 31) * self._stoichs['27_10']
        self._stoichs['27_22'] = self._params['i_TSSBM'] * self._stoichs['27_17'] + 0.6 * self._stoichs['27_19']

        self._stoichs['28_4'] = self._params['i_NBM'] - self._params['i_NXI'] * self._params['f_XI'] - (1.0 - self._params['f_XI']) * self._params['i_NXS']
        self._stoichs['28_10'] = self._params['i_PBM'] - self._params['i_PXI'] * self._params['f_XI'] - (1.0 - self._params['f_XI']) * self._params['i_PXS']
        self._stoichs['28_14'] = self._params['f_XI']
        self._stoichs['28_15'] = 1.0 - self._params['f_XI']
        self._stoichs['28_17'] = -1.0
        self._stoichs['28_12'] = (1 / 14) * self._stoichs['28_4'] + (-1.5 / 31) * self._stoichs['28_10']
        self._stoichs['28_22'] = self._params['i_TSSXI'] * self._stoichs['28_14'] + self._params['i_TSSXS'] * self._stoichs['28_15'] + self._params['i_TSSBM'] * self._stoichs['28_17']

        self._stoichs['29_10'] = 1.0
        self._stoichs['29_18'] = -1.0
        self._stoichs['29_12'] = (-1.5 / 31) * self._stoichs['29_10'] + (-1 / 31) * self._stoichs['29_18']
        self._stoichs['29_22'] = 3.23 * self._stoichs['29_18']

        self._stoichs['30_3'] = 1.0
        self._stoichs['30_19'] = -1.0
        self._stoichs['30_12'] = (-1 / 64) * self._stoichs['30_3']
        self._stoichs['30_22'] = 0.6 * self._stoichs['30_19']

        self._stoichs['31_1'] = -8 / 7
        self._stoichs['31_4'] = -1.0
        self._stoichs['31_5'] = 1.0
        self._stoichs['31_12'] = (1 / 14) * self._stoichs['31_4']

        self._stoichs['32_1'] = 0.0 - ((12 / 7) - self._params['Y_AOB']) / self._params['Y_AOB']
        self._stoichs['32_4'] = 0.0 - self._params['i_NBM']
        self._stoichs['32_5'] = -1.0 / self._params['Y_AOB']
        self._stoichs['32_7'] = 1.0 / self._params['Y_AOB']
        self._stoichs['32_10'] = 0.0 - self._params['i_PBM']
        self._stoichs['32_20'] = 1.0
        self._stoichs['32_12'] = (1 / 14) * self._stoichs['32_4'] + (-1.5 / 31) * self._stoichs['32_10']
        self._stoichs['32_22'] = self._params['i_TSSBM'] * self._stoichs['32_20']

        self._stoichs['33_1'] = -4 / 7
        self._stoichs['33_7'] = -1.0
        self._stoichs['33_8'] = 1.0
        self._stoichs['33_12'] = (-1 / 14) * self._stoichs['33_8']

        self._stoichs['34_5'] = -1.0
        self._stoichs['34_6'] = 4.0
        self._stoichs['34_7'] = -4.0
        self._stoichs['34_8'] = 1.0
        self._stoichs['34_12'] = (-1 / 14) * self._stoichs['34_8']

        self._stoichs['35_5'] = -1.0
        self._stoichs['35_6'] = 2.0
        self._stoichs['35_8'] = -1.0
        self._stoichs['35_12'] = (-1 / 14) * self._stoichs['35_8']

        self._stoichs['36_1'] = 0.0 - ((8 / 7) - self._params['Y_NOB']) / self._params['Y_NOB']
        self._stoichs['36_4'] = 0.0 - self._params['i_NBM']
        self._stoichs['36_8'] = -1.0 / self._params['Y_NOB']
        self._stoichs['36_9'] = 1.0 / self._params['Y_NOB']
        self._stoichs['36_10'] = 0.0 - self._params['i_PBM']
        self._stoichs['36_21'] = 1.0
        self._stoichs['36_12'] = (1 / 14) * self._stoichs['36_4'] + (-1 / 14) * self._stoichs['36_8'] + (-1 / 14) * self._stoichs['36_9'] + (-1.5 / 31) * self._stoichs['36_10']
        self._stoichs['36_22'] = self._params['i_TSSBM'] * self._stoichs['36_21']

        self._stoichs['37_4'] = self._params['i_NBM'] - self._params['i_NXI'] * self._params['f_XI'] - (1.0 - self._params['f_XI']) * self._params['i_NXS']
        self._stoichs['37_10'] = self._params['i_PBM'] - self._params['i_PXI'] * self._params['f_XI'] - (1.0 - self._params['f_XI']) * self._params['i_PXS']
        self._stoichs['37_14'] = self._params['f_XI']
        self._stoichs['37_15'] = 1.0 - self._params['f_XI']
        self._stoichs['37_20'] = -1.0
        self._stoichs['37_12'] = (1 / 14) * self._stoichs['37_4'] + (-1.5 / 31) * self._stoichs['37_10']
        self._stoichs['37_22'] = self._params['i_TSSXI'] * self._stoichs['37_14'] + self._params['i_TSSXS'] * self._stoichs['37_15'] + self._params['i_TSSBM'] * self._stoichs['37_20']
        
        self._stoichs['38_4'] = self._params['i_NBM'] - self._params['i_NXI'] * self._params['f_XI'] - (1.0 - self._params['f_XI']) * self._params['i_NXS']
        self._stoichs['38_10'] = self._params['i_PBM'] - self._params['i_PXI'] * self._params['f_XI'] - (1.0 - self._params['f_XI']) * self._params['i_PXS']
        self._stoichs['38_14'] = self._params['f_XI']
        self._stoichs['38_15'] = 1.0 - self._params['f_XI']
        self._stoichs['38_21'] = -1.0
        self._stoichs['38_12'] = (1 / 14) * self._stoichs['38_4'] + (-1.5 / 31) * self._stoichs['38_10']
        self._stoichs['38_22'] = self._params['i_TSSXI'] * self._stoichs['38_14'] + self._params['i_TSSXS'] * self._stoichs['38_15'] + self._params['i_TSSBM'] * self._stoichs['38_21']

        self._stoichs['39_10'] = -1.0
        self._stoichs['39_23'] = -3.45
        self._stoichs['39_24'] = 4.87
        self._stoichs['39_12'] = (-1.5 / 31) * self._stoichs['39_10']
        self._stoichs['39_22'] = 1.0 * self._stoichs['39_23'] + 1.0 * self._stoichs['39_24']

        self._stoichs['40_10'] = 1.0
        self._stoichs['40_23'] = 3.45
        self._stoichs['40_24'] = -4.87
        self._stoichs['40_12'] = (-1.5 / 31) * self._stoichs['40_10']
        self._stoichs['40_22'] = 1.0 * self._stoichs['40_23'] + 1.0 * self._stoichs['40_24']

        return None
    
    def _reaction_rate(self, comps):
        '''
        Calculate the reaction rate for each biological process.

        Args:
            comps (list): A list of current model components (concentrations).

        Return:
            self._rate_res (list): A list of reaction rates for each biological process.
        '''

        ## Use some Monod terms for simplification
        self._monods[0] = self._monod(comps[0], self._params['K_O2_H'])
        self._monods[1] = self._monod((comps[14] / comps[15]), self._params['K_X_H'])
        self._monods[2] = self._monod(self._params['K_O2_H'], comps[0])
        self._monods[3] = self._monod(comps[8], self._params['K_NO3_H'])
        self._monods[4] = self._monod(comps[7], self._params['K_NO2_H'])
        self._monods[5] = self._monod(self._params['K_NO2_H'], (comps[7] + comps[8]))
        self._monods[6] = self._monod(comps[1], self._params['K_F'])
        self._monods[7] = self._monod(comps[1], comps[2])
        self._monods[8] = self._monod(comps[0], self._params['K_O2'])
        self._monods[9] = self._monod(comps[3], self._params['K_NH4'])
        self._monods[10] = self._monod(comps[9], self._params['K_P'])
        self._monods[11] = self._monod(comps[11], self._params['K_ALK'])
        self._monods[12] = self._monod(comps[2], self._params['K_A'])
        self._monods[13] = self._monod(comps[2], comps[1])
        self._monods[14] = self._monod(self._params['K_O2'], comps[0])
        self._monods[15] = self._monod(comps[1], self._params['K_S3'])
        self._monods[16] = self._monod(comps[7], self._params['K_NO2_Den'])
        self._monods[17] = self._monod(self._params['K_OH3'], comps[0])
        self._monods[18] = self._monod(comps[1], self._params['K_S4'])
        self._monods[19] = self._monod(comps[6], (self._params['K_NO_Den'] + (comps[6])**2 / self._params['K_I4NO']))
        self._monods[20] = self._monod(self._params['K_OH4'], comps[0])
        self._monods[21] = self._monod(comps[1], self._params['K_S5'])
        self._monods[22] = self._monod(comps[5], self._params['K_N2O_Den'])
        self._monods[23] = self._monod(self._params['K_OH5'], comps[0])
        self._monods[24] = self._monod(comps[2], self._params['K_S3'])
        self._monods[25] = self._monod(comps[2], self._params['K_S4'])
        self._monods[26] = self._monod(comps[2], self._params['K_S5'])
        self._monods[27] = self._monod(self._params['K_NO2'], (comps[7] + comps[8]))
        self._monods[28] = self._monod(comps[1], self._params['K_fe_H'])
        self._monods[29] = self._monod(comps[2], self._params['K_A_P'])
        self._monods[30] = self._monod(comps[11], self._params['K_ALK_P'])
        self._monods[31] = self._monod((comps[17] / comps[16]), self._params['K_PP_P'])
        self._monods[32] = self._monod(comps[0], self._params['K_O2_P'])
        self._monods[33] = self._monod(comps[9], self._params['K_P_P'])
        self._monods[34] = self._monod((comps[18] / comps[16]), self._params['K_PHA_P'])
        self._monods[35] = self._monod((self._params['K_MAX_P'] - (comps[17] / comps[16])), self._params['K_IPP_P'])
        self._monods[36] = self._monod(comps[8], self._params['K_NO3_P'])
        self._monods[37] = self._monod(self._params['K_O2_P'], comps[0])
        self._monods[38] = self._monod(comps[0], self._params['K_O2_AOB1'])
        self._monods[39] = self._monod(comps[3], self._params['K_NH4_AOB'])
        self._monods[40] = self._monod(comps[0], self._params['K_O2_AOB2'])
        self._monods[41] = self._monod(comps[4], self._params['K_NH2OH_AOB'])
        self._monods[42] = self._monod(comps[9], self._params['K_P_AOB'])
        self._monods[43] = self._monod(comps[11], self._params['K_ALK_AOB'])
        self._monods[44] = self._monod(comps[6], self._params['K_NO_AOB_HAO'])
        self._monods[45] = self._monod(comps[6], self._params['K_NO_AOB_NN'])
        
        # Calculate concentrations of HNO2
        temp = 20
        pH = 7
        Ka = math.exp(-2300 / (273.15 + temp))
        S_HNO2 = (comps[7] / (Ka * 10**pH + 1)) * (47 / 14)
        self._monods[46] = self._monod(S_HNO2, self._params['K_HNO2_AOB'])

        self._monods[47] = self._monod(comps[0], self._params['K_O2_NOB'])
        self._monods[48] = self._monod(comps[7], self._params['K_NO2_NOB'])
        self._monods[49] = self._monod(comps[9], self._params['K_P_NOB'])
        self._monods[50] = self._monod(comps[11], self._params['K_ALK_NOB'])
        self._monods[51] = self._monod(comps[11], self._params['K_ALK_PR'])

        # Omitted monod terms earlier
        self._monods[52] = self._monod(comps[8], self._params['K_NO3'])
        self._monods[53] = self._monod(comps[3], 10e-12)

        ## Calculate reaction rates
        self._rate_res[0] = self._params['K_H'] * self._monods[0] * self._monods[1] * comps[15]
        self._rate_res[1] = self._params['K_H'] * self._params['n_NO3_H'] * self._monods[2] * self._monods[3] * self._monods[1] * comps[15]
        self._rate_res[2] = self._params['K_H'] * self._params['n_NO2_H'] * self._monods[2] * self._monods[4] * self._monods[1] * comps[15]
        self._rate_res[3] = self._params['K_H'] * self._params['n_fe_H'] * self._monods[2] * self._monods[5] * self._monods[1] * comps[15]
        self._rate_res[4] = self._params['mu_H'] * self._monods[6] * self._monods[7] * self._monods[8] * self._monods[9] * self._monods[10] * self._monods[11] * comps[15]
        self._rate_res[5] = self._params['mu_H'] * self._monods[12] * self._monods[13] * self._monods[8] * self._monods[9] * self._monods[10] * self._monods[11] * comps[15]
        self._rate_res[6] = self._params['mu_H'] * self._params['n_NO3_D'] * self._monods[6] * self._monods[7] * self._monods[14] * self._monods[52] * self._monods[9] * self._monods[10] * self._monods[11] * comps[15]
        self._rate_res[7] = self._params['mu_H'] * self._params['n_G3'] * self._monods[15] * self._monods[7] * self._monods[16] * self._monods[17] * self._monods[9] * self._monods[10] * self._monods[11] * comps[15]
        self._rate_res[8] = self._params['mu_H'] * self._params['n_G4'] * self._monods[18] * self._monods[7] * self._monods[19] * self._monods[20] * self._monods[9] * self._monods[10] * self._monods[11] * comps[15]
        self._rate_res[9] = self._params['mu_H'] * self._params['n_G5'] * self._monods[21] * self._monods[7] * self._monods[22] * self._monods[23] * self._monods[9] * self._monods[10] * self._monods[11] * comps[15]
        self._rate_res[10] = self._params['mu_H'] * self._params['n_NO3_D'] * self._monods[12] * self._monods[13] * self._monods[14] * self._monods[52] * self._monods[9] * self._monods[10] * self._monods[11] * comps[15]
        self._rate_res[11] = self._params['mu_H'] * self._params['n_G3'] * self._monods[24] * self._monods[13] * self._monods[16] * self._monods[17] * self._monods[9] * self._monods[10] * self._monods[11] * comps[15]
        self._rate_res[12] = self._params['mu_H'] * self._params['n_G4'] * self._monods[25] * self._monods[13] * self._monods[19] * self._monods[20] * self._monods[9] * self._monods[10] * self._monods[11] * comps[15]
        self._rate_res[13] = self._params['mu_H'] * self._params['n_G5'] * self._monods[26] * self._monods[13] * self._monods[22] * self._monods[23] * self._monods[9] * self._monods[10] * self._monods[11] * comps[15]
        self._rate_res[14] = self._params['q_fe'] * self._monods[14] * self._monods[27] * self._monod[28] * self._monods[11] * comps[15]
        self._rate_res[15] = self._params['b_H'] * comps[15]
        self._rate_res[16] = self._params['q_PHA'] * self._monods[29] * self._monods[30] * self._monods[31] * comps[16]
        self._rate_res[17] = self._params['q_PP'] * self._monods[32] * self._monods[33] * self._monods[30] * self._monods[34] * self._monods[35] * comps[16]
        self._rate_res[18] = self._params['q_PP'] * self._params['n_NO3_P'] * self._monods[36] * self._monods[37] * self._monods[33] * self._monods[30] * self._monods[34] * self._monods[35] * comps[16]
        self._rate_res[19] = self._params['q_PP'] * self._params['n_G3'] * self._monods[16] * self._monods[17] * self._monods[33] * self._monods[30] * self._monods[34] * self._monods[35] * comps[16]
        self._rate_res[20] = self._params['q_PP'] * self._params['n_G4'] * self._monods[19] * self._monods[20] * self._monods[33] * self._monods[30] * self._monods[34] * self._monods[35] * comps[16]
        self._rate_res[21] = self._params['q_PP'] * self._params['n_G5'] * self._monods[22] * self._monods[23] * self._monods[33] * self._monods[30] * self._monods[34] * self._monods[35] * comps[16]
        self._rate_res[22] = self._params['mu_PAO'] * self._monods[32] * self._monods[33] * self._monods[9] * self._monods[30] * self._monods[34] * comps[16]
        self._rate_res[23] = self._params['mu_PAO'] * self._params['n_NO3_P'] * self._monods[36] * self._monods[37] * self._monods[33] * self._monods[9] * self._monods[30] * self._monods[34] * comps[16]
        self._rate_res[24] = self._params['mu_PAO'] * self._params['n_G3'] * self._monods[16] * self._monods[17] * self._monods[33] * self._monods[9] * self._monods[30] * self._monods[34] * comps[16]
        self._rate_res[25] = self._params['mu_PAO'] * self._params['n_G4'] * self._monods[19] * self._monods[20] * self._monods[33] * self._monods[9] * self._monods[30] * self._monods[34] * comps[16]
        self._rate_res[26] = self._params['mu_PAO'] * self._params['n_G5'] * self._monods[22] * self._monods[23] * self._monods[33] * self._monods[9] * self._monods[30] * self._monods[34] * comps[16]
        self._rate_res[27] = self._params['b_PAO'] * self._monods[30] * comps[16]
        self._rate_res[28] = self._params['b_PP'] * self._monods[30] * comps[17]
        self._rate_res[29] = self._params['b_PHA'] * self._monods[30] * comps[18]
        self._rate_res[30] = self._params['q_AOB_AMO'] * self._monods[38] * self._monods[39] * comps[19]
        self._rate_res[31] = self._params['mu_AOB_HAO'] * self._monods[40] * self._monods[41] * self._monods[53] * self._monods[42] * self._monods[43] * comps[19]
        self._rate_res[32] = self._params['q_AOB_HAO'] * self._monods[40] * self._monods[44] * comps[19]
        self._rate_res[33] = self._params['q_AOB_N2O_NN'] * self._monods[41] * self._monods[45] * comps[19]

        # Calculate concentrations of HNO2
        fSO2 = comps[0] / (self._params['K_O2_AOB_ND'] + (1.0 - 2.0 * (self._params['K_O2_AOB_ND'] / self._params['K_I_O2_AOB']) ** (1 / 2)) * comps[0] + ((comps[0] ** 2) / self._params['K_I_O2_AOB']))
        self._rate_res[34] = self._params['q_AOB_N2O_ND'] * self._monods[41] * self._monods[46] * fSO2 * comps[19]

        self._rate_res[35] = self._params['mu_NOB'] * self._monods[47] * self._monods[48] * self._monods[49] * self._monods[50] * comps[20]
        self._rate_res[36] = self._params['b_AOB'] * comps[19]
        self._rate_res[37] = self._params['b_NOB'] * comps[20]
        self._rate_res[38] = self._params['k_PRE'] * comps[9] * comps[22]
        self._rate_res[39] = self._params['k_RED'] * self._monods[51] * comps[23]

        return self._rate_res
    
    ## Overall process rates for each component

    def _rate0_S_O2(self):
        '''
        Overall process rate for dissolved oxygen (mgCOD/L/d).
        
        Return:
            rate (float): Process rate for dissolved oxygen.
        '''
        return (self._stoichs['5_1'] * self._rate_res[4]
                + self._stoichs['6_1'] * self._rate_res[5]
                + self._stoichs['18_1'] * self._rate_res[17]
                + self._stoichs['23_1'] * self._rate_res[22]
                + self._stoichs['31_1'] * self._rate_res[30]
                + self._stoichs['32_1'] * self._rate_res[31]
                + self._stoichs['33_1'] * self._rate_res[32]
                + self._stoichs['36_1'] * self._rate_res[35])

    def _rate1_S_F(self):
        '''
        Overall process rate for fermentable, readily biodegradable organic substrates (mgCOD/L/d).
        
        Return:
            rate (float): Process rate for fermentable, readily biodegradable organic substrates.
        '''
        return (self._stoichs['1_2'] * self._rate_res[0]
                + self._stoichs['2_2'] * self._rate_res[1]
                + self._stoichs['3_2'] * self._rate_res[2]
                + self._stoichs['4_2'] * self._rate_res[3]
                + self._stoichs['5_2'] * self._rate_res[4]
                + self._stoichs['7_2'] * self._rate_res[6]
                + self._stoichs['8_2'] * self._rate_res[7]
                + self._stoichs['9_2'] * self._rate_res[8]
                + self._stoichs['10_2'] * self._rate_res[9]
                + self._stoichs['15_2'] * self._rate_res[14])

    def _rate2_S_A(self):
        '''
        Overall process rate for acetate (mgCOD/L/d).
        
        Return:
            rate (float): Process rate for acetate.
        '''
        return (self._stoichs['6_3'] * self._rate_res[5]
                + self._stoichs['11_3'] * self._rate_res[10]
                + self._stoichs['12_3'] * self._rate_res[11]
                + self._stoichs['13_3'] * self._rate_res[12]
                + self._stoichs['14_3'] * self._rate_res[13]
                + self._stoichs['15_3'] * self._rate_res[14]
                + self._stoichs['17_3'] * self._rate_res[16]
                + self._stoichs['30_3'] * self._rate_res[29])

    def _rate3_S_NH4(self):
        '''
        Overall process rate for ammonium (mgN/L/d).
        
        Return:
            rate (float): Process rate for ammonium.
        '''
        return (self._stoichs['1_4'] * self._rate_res[0]
                + self._stoichs['2_4'] * self._rate_res[1]
                + self._stoichs['3_4'] * self._rate_res[2]
                + self._stoichs['4_4'] * self._rate_res[3]
                + self._stoichs['5_4'] * self._rate_res[4]
                + self._stoichs['6_4'] * self._rate_res[5]
                + self._stoichs['7_4'] * self._rate_res[6]
                + self._stoichs['8_4'] * self._rate_res[7]
                + self._stoichs['9_4'] * self._rate_res[8]
                + self._stoichs['10_4'] * self._rate_res[9]
                + self._stoichs['11_4'] * self._rate_res[10]
                + self._stoichs['12_4'] * self._rate_res[11]
                + self._stoichs['13_4'] * self._rate_res[12]
                + self._stoichs['14_4'] * self._rate_res[13]
                + self._stoichs['15_4'] * self._rate_res[14]
                + self._stoichs['16_4'] * self._rate_res[15]
                + self._stoichs['23_4'] * self._rate_res[22]
                + self._stoichs['24_4'] * self._rate_res[23]
                + self._stoichs['25_4'] * self._rate_res[24]
                + self._stoichs['26_4'] * self._rate_res[25]
                + self._stoichs['27_4'] * self._rate_res[26]
                + self._stoichs['28_4'] * self._rate_res[27]
                + self._stoichs['31_4'] * self._rate_res[30]
                + self._stoichs['32_4'] * self._rate_res[31]
                + self._stoichs['36_4'] * self._rate_res[35]
                + self._stoichs['37_4'] * self._rate_res[36]
                + self._stoichs['38_4'] * self._rate_res[37])

    def _rate4_S_NH2OH(self):
        '''
        Overall process rate for hydroxylamine (mgN/L/d).
        
        Return:
            rate (float): Process rate for hydroxylamine.
        '''
        return (self._stoichs['31_5'] * self._rate_res[30]
                + self._stoichs['32_5'] * self._rate_res[31]
                + self._stoichs['34_5'] * self._rate_res[33]
                + self._stoichs['35_5'] * self._rate_res[34])
    
    def _rate5_S_N2O(self):
        '''
        Overall process rate for nitrous oxide (mgN/L/d).
        
        Return:
            rate (float): Process rate for nitrous oxide.
        '''
        return (self._stoichs['9_6'] * self._rate_res[8]
                + self._stoichs['10_6'] * self._rate_res[9]
                + self._stoichs['13_6'] * self._rate_res[12]
                + self._stoichs['14_6'] * self._rate_res[13]
                + self._stoichs['21_6'] * self._rate_res[20]
                + self._stoichs['22_6'] * self._rate_res[21]
                + self._stoichs['26_6'] * self._rate_res[25]
                + self._stoichs['27_6'] * self._rate_res[26]
                + self._stoichs['34_6'] * self._rate_res[33]
                + self._stoichs['35_6'] * self._rate_res[34])
    
    def _rate6_S_NO(self):
        '''
        Overall process rate for nitric oxide (mgN/L/d).
        
        Return:
            rate (float): Process rate for nitric oxide.
        '''
        return (self._stoichs['8_7'] * self._rate_res[7]
                + self._stoichs['9_7'] * self._rate_res[8]
                + self._stoichs['12_7'] * self._rate_res[11]
                + self._stoichs['13_7'] * self._rate_res[12]
                + self._stoichs['20_7'] * self._rate_res[19]
                + self._stoichs['21_7'] * self._rate_res[20]
                + self._stoichs['25_7'] * self._rate_res[24]
                + self._stoichs['26_7'] * self._rate_res[25]
                + self._stoichs['32_7'] * self._rate_res[31]
                + self._stoichs['33_7'] * self._rate_res[32]
                + self._stoichs['34_7'] * self._rate_res[33])
    
    def _rate7_S_NO2(self):
        '''
        Overall process rate for nitrite (mgN/L/d).
        
        Return:
            rate (float): Process rate for nitrite.
        '''
        return (self._stoichs['7_8'] * self._rate_res[6]
                + self._stoichs['8_8'] * self._rate_res[7]
                + self._stoichs['11_8'] * self._rate_res[10]
                + self._stoichs['12_8'] * self._rate_res[11]
                + self._stoichs['19_8'] * self._rate_res[18]
                + self._stoichs['20_8'] * self._rate_res[19]
                + self._stoichs['24_8'] * self._rate_res[23]
                + self._stoichs['25_8'] * self._rate_res[24]
                + self._stoichs['33_8'] * self._rate_res[32]
                + self._stoichs['34_8'] * self._rate_res[33]
                + self._stoichs['35_8'] * self._rate_res[34]
                + self._stoichs['36_8'] * self._rate_res[35])
    
    def _rate8_S_NO3(self):
        '''
        Overall process rate for nitrate (mgN/L/d).
        
        Return:
            rate (float): Process rate for nitrate.
        '''
        return (self._stoichs['7_9'] * self._rate_res[6]
                + self._stoichs['11_9'] * self._rate_res[10]
                + self._stoichs['19_9'] * self._rate_res[18]
                + self._stoichs['24_9'] * self._rate_res[23]
                + self._stoichs['36_9'] * self._rate_res[35])
    
    def _rate9_S_PO4(self):
        '''
        Overall process rate for phosphate (mgP/L/d).
        
        Return:
            rate (float): Process rate for phosphate.
        '''
        return (self._stoichs['1_10'] * self._rate_res[0]
                + self._stoichs['2_10'] * self._rate_res[1]
                + self._stoichs['3_10'] * self._rate_res[2]
                + self._stoichs['4_10'] * self._rate_res[3]
                + self._stoichs['5_10'] * self._rate_res[4]
                + self._stoichs['6_10'] * self._rate_res[5]
                + self._stoichs['7_10'] * self._rate_res[6]
                + self._stoichs['8_10'] * self._rate_res[7]
                + self._stoichs['9_10'] * self._rate_res[8]
                + self._stoichs['10_10'] * self._rate_res[9]
                + self._stoichs['11_10'] * self._rate_res[10]
                + self._stoichs['12_10'] * self._rate_res[11]
                + self._stoichs['13_10'] * self._rate_res[12]
                + self._stoichs['14_10'] * self._rate_res[13]
                + self._stoichs['15_10'] * self._rate_res[14]
                + self._stoichs['16_10'] * self._rate_res[15]
                + self._stoichs['17_10'] * self._rate_res[16]
                + self._stoichs['18_10'] * self._rate_res[17]
                + self._stoichs['19_10'] * self._rate_res[18]
                + self._stoichs['20_10'] * self._rate_res[19]
                + self._stoichs['21_10'] * self._rate_res[20]
                + self._stoichs['22_10'] * self._rate_res[21]
                + self._stoichs['23_10'] * self._rate_res[22]
                + self._stoichs['24_10'] * self._rate_res[23]
                + self._stoichs['25_10'] * self._rate_res[24]
                + self._stoichs['26_10'] * self._rate_res[25]
                + self._stoichs['27_10'] * self._rate_res[26]
                + self._stoichs['28_10'] * self._rate_res[27]
                + self._stoichs['29_10'] * self._rate_res[28]
                + self._stoichs['32_10'] * self._rate_res[31]
                + self._stoichs['36_10'] * self._rate_res[35]
                + self._stoichs['37_10'] * self._rate_res[36]
                + self._stoichs['38_10'] * self._rate_res[37]
                + self._stoichs['39_10'] * self._rate_res[38]
                + self._stoichs['40_10'] * self._rate_res[39])
    
    def _rate10_S_I(self):
        '''
        Overall process rate for inert soluble material (mgCOD/L/d).
        
        Return:
            rate (float): Process rate for inert soluble material.
        '''
        return (self._stoichs['1_11'] * self._rate_res[0]
                + self._stoichs['2_11'] * self._rate_res[1]
                + self._stoichs['3_11'] * self._rate_res[2]
                + self._stoichs['4_11'] * self._rate_res[3])
    
    def _rate11_S_ALK(self):
        '''
        Overall process rate for alkalinity (molHCO3/m3/d).
        
        Return:
            rate (float): Process rate for alkalinity.
        '''
        return (self._stoichs['1_12'] * self._rate_res[0]
                + self._stoichs['2_12'] * self._rate_res[1]
                + self._stoichs['3_12'] * self._rate_res[2]
                + self._stoichs['4_12'] * self._rate_res[3]
                + self._stoichs['5_12'] * self._rate_res[4]
                + self._stoichs['6_12'] * self._rate_res[5]
                + self._stoichs['7_12'] * self._rate_res[6]
                + self._stoichs['8_12'] * self._rate_res[7]
                + self._stoichs['9_12'] * self._rate_res[8]
                + self._stoichs['10_12'] * self._rate_res[9]
                + self._stoichs['11_12'] * self._rate_res[10]
                + self._stoichs['12_12'] * self._rate_res[11]
                + self._stoichs['13_12'] * self._rate_res[12]
                + self._stoichs['14_12'] * self._rate_res[13]
                + self._stoichs['15_12'] * self._rate_res[14]
                + self._stoichs['16_12'] * self._rate_res[15]
                + self._stoichs['17_12'] * self._rate_res[16]
                + self._stoichs['18_12'] * self._rate_res[17]
                + self._stoichs['19_12'] * self._rate_res[18]
                + self._stoichs['20_12'] * self._rate_res[19]
                + self._stoichs['21_12'] * self._rate_res[20]
                + self._stoichs['22_12'] * self._rate_res[21]
                + self._stoichs['23_12'] * self._rate_res[22]
                + self._stoichs['24_12'] * self._rate_res[23]
                + self._stoichs['25_12'] * self._rate_res[24]
                + self._stoichs['26_12'] * self._rate_res[25]
                + self._stoichs['27_12'] * self._rate_res[26]
                + self._stoichs['28_12'] * self._rate_res[27]
                + self._stoichs['29_12'] * self._rate_res[28]
                + self._stoichs['30_12'] * self._rate_res[29]
                + self._stoichs['31_12'] * self._rate_res[30]
                + self._stoichs['32_12'] * self._rate_res[31]
                + self._stoichs['33_12'] * self._rate_res[32]
                + self._stoichs['34_12'] * self._rate_res[33]
                + self._stoichs['35_12'] * self._rate_res[34]
                + self._stoichs['36_12'] * self._rate_res[35]
                + self._stoichs['37_12'] * self._rate_res[36]
                + self._stoichs['38_12'] * self._rate_res[37]
                + self._stoichs['39_12'] * self._rate_res[38]
                + self._stoichs['40_12'] * self._rate_res[39])
    
    def _rate12_S_N2(self):
        '''
        Overall process rate for nitrogen gas (mgN/L/d).
        
        Return:
            rate (float): Process rate for nitrogen gas.
        '''
        return (self._stoichs['10_13'] * self._rate_res[9]
                + self._stoichs['14_13'] * self._rate_res[13]
                + self._stoichs['22_13'] * self._rate_res[21]
                + self._stoichs['27_13'] * self._rate_res[26])

    def _rate13_X_I(self):
        '''
        Overall process rate for inert particulate organic material (mgCOD/L/d).
        
        Return:
            rate (float): Process rate for inert particulate organic material.
        '''
        return (self._stoichs['16_14'] * self._rate_res[15]
                + self._stoichs['28_14'] * self._rate_res[27]
                + self._stoichs['37_14'] * self._rate_res[36]
                + self._stoichs['38_14'] * self._rate_res[37])
    
    def _rate14_X_S(self):
        '''
        Overall process rate for slowly biodegradable substrates (mgCOD/L/d).
        
        Return:
            rate (float): Process rate for slowly biodegradable substrates.
        '''
        return (self._stoichs['1_15'] * self._rate_res[0]
                + self._stoichs['2_15'] * self._rate_res[1]
                + self._stoichs['3_15'] * self._rate_res[2]
                + self._stoichs['4_15'] * self._rate_res[3]
                + self._stoichs['16_15'] * self._rate_res[15]
                + self._stoichs['28_15'] * self._rate_res[27]
                + self._stoichs['37_15'] * self._rate_res[36]
                + self._stoichs['38_15'] * self._rate_res[37])
    
    def _rate15_X_H(self):
        '''
        Overall process rate for heterotrophic biomass (mgCOD/L/d).
        
        Return:
            rate (float): Process rate for heterotrophic biomass.
        '''
        return (self._stoichs['5_16'] * self._rate_res[4]
                + self._stoichs['6_16'] * self._rate_res[5]
                + self._stoichs['7_16'] * self._rate_res[6]
                + self._stoichs['8_16'] * self._rate_res[7]
                + self._stoichs['9_16'] * self._rate_res[8]
                + self._stoichs['10_16'] * self._rate_res[9]
                + self._stoichs['11_16'] * self._rate_res[10]
                + self._stoichs['12_16'] * self._rate_res[11]
                + self._stoichs['13_16'] * self._rate_res[12]
                + self._stoichs['14_16'] * self._rate_res[13]
                + self._stoichs['16_16'] * self._rate_res[15])
    
    def _rate16_X_PAO(self):
        '''
        Overall process rate for P-accumulating organisms (mgCOD/L/d).
        
        Return:
            rate (float): Process rate for P-accumulating organisms.
        '''
        return (self._stoichs['23_17'] * self._rate_res[22]
                + self._stoichs['24_17'] * self._rate_res[23]
                + self._stoichs['25_17'] * self._rate_res[24]
                + self._stoichs['26_17'] * self._rate_res[25]
                + self._stoichs['27_17'] * self._rate_res[26]
                + self._stoichs['28_17'] * self._rate_res[27])
    
    def _rate17_X_PP(self):
        '''
        Overall process rate for Polyphosphate Accumulating Organisms (mgCOD/L/d).
        
        Return:
            rate (float): Process rate for Polyphosphate Accumulating Organisms.
        '''
        return (self._stoichs['17_18'] * self._rate_res[16]
                + self._stoichs['18_18'] * self._rate_res[17]
                + self._stoichs['19_18'] * self._rate_res[18]
                + self._stoichs['20_18'] * self._rate_res[19]
                + self._stoichs['21_18'] * self._rate_res[20]
                + self._stoichs['22_18'] * self._rate_res[21]
                + self._stoichs['29_18'] * self._rate_res[28])

    def _rate18_X_PHA(self):
        '''
        Overall process rate for polyhydroxyalkanoates accumulating organisms (mgCOD/L/d).
        
        Return:
            rate (float): Process rate for polyhydroxyalkanoates accumulating organisms.
        '''
        return (self._stoichs['17_19'] * self._rate_res[16]
                + self._stoichs['18_19'] * self._rate_res[17]
                + self._stoichs['19_19'] * self._rate_res[18]
                + self._stoichs['20_19'] * self._rate_res[19]
                + self._stoichs['21_19'] * self._rate_res[20]
                + self._stoichs['22_19'] * self._rate_res[21]
                + self._stoichs['23_19'] * self._rate_res[22]
                + self._stoichs['24_19'] * self._rate_res[23]
                + self._stoichs['25_19'] * self._rate_res[24]
                + self._stoichs['26_19'] * self._rate_res[25]
                + self._stoichs['27_19'] * self._rate_res[26]
                + self._stoichs['30_19'] * self._rate_res[29])
    
    def _rate19_X_AOB(self):
        '''
        Overall process rate for ammonia oxidizing bacteria (mgCOD/L/d).
        
        Return:
            rate (float): Process rate for ammonia oxidizing bacteria.
        '''
        return (self._stoichs['32_20'] * self._rate_res[31]
                + self._stoichs['37_20'] * self._rate_res[36])
    
    def _rate20_X_NOB(self):
        '''
        Overall process rate for nitrite oxidizing bacteria (mgCOD/L/d).
        
        Return:
            rate (float): Process rate for nitrite oxidizing bacteria.
        '''
        return (self._stoichs['36_21'] * self._rate_res[35]
                + self._stoichs['38_21'] * self._rate_res[37])

    def _rate21_X_TSS(self):
        '''
        Overall process rate for total suspended solids (mgCOD/L/d).

        Return:
            rate (float): Process rate for total suspended solids.
        '''
        return (self._stoichs['1_22'] * self._rate_res[0]
                + self._stoichs['2_22'] * self._rate_res[1]
                + self._stoichs['3_22'] * self._rate_res[2]
                + self._stoichs['4_22'] * self._rate_res[3]
                + self._stoichs['5_22'] * self._rate_res[4]
                + self._stoichs['6_22'] * self._rate_res[5]
                + self._stoichs['7_22'] * self._rate_res[6]
                + self._stoichs['8_22'] * self._rate_res[7]
                + self._stoichs['9_22'] * self._rate_res[8]
                + self._stoichs['10_22'] * self._rate_res[9]
                + self._stoichs['11_22'] * self._rate_res[10]
                + self._stoichs['12_22'] * self._rate_res[11]
                + self._stoichs['13_22'] * self._rate_res[12]
                + self._stoichs['14_22'] * self._rate_res[13]
                + self._stoichs['16_22'] * self._rate_res[15]
                + self._stoichs['17_22'] * self._rate_res[16]
                + self._stoichs['18_22'] * self._rate_res[17]
                + self._stoichs['19_22'] * self._rate_res[18]
                + self._stoichs['20_22'] * self._rate_res[19]
                + self._stoichs['21_22'] * self._rate_res[20]
                + self._stoichs['22_22'] * self._rate_res[21]
                + self._stoichs['23_22'] * self._rate_res[22]
                + self._stoichs['24_22'] * self._rate_res[23]
                + self._stoichs['25_22'] * self._rate_res[24]
                + self._stoichs['26_22'] * self._rate_res[25]
                + self._stoichs['27_22'] * self._rate_res[26]
                + self._stoichs['28_22'] * self._rate_res[27]
                + self._stoichs['29_22'] * self._rate_res[28]
                + self._stoichs['30_22'] * self._rate_res[29]
                + self._stoichs['32_22'] * self._rate_res[31]
                + self._stoichs['36_22'] * self._rate_res[35]
                + self._stoichs['37_22'] * self._rate_res[36]
                + self._stoichs['38_22'] * self._rate_res[37]
                + self._stoichs['39_22'] * self._rate_res[38]
                + self._stoichs['40_22'] * self._rate_res[39])

    def _rate22_X_MeOH(self):
        '''
        Overall process rate for P-precipitation with Fe(OH)3 (mgP/L/d).
        
        Return:
            rate (float): Process rate for P-precipitation.
        '''
        return (self._stoichs['39_23'] * self._rate_res[38]
                + self._stoichs['40_23'] * self._rate_res[39])
    
    def _rate23_X_MeP(self):
        '''
        Overall process rate for P-redissolution (mgP/L/d).
        
        Return:
            rate (float): Process rate for P-redissolution.
        '''
        return (self._stoichs['39_24'] * self._rate_res[38]
                + self._stoichs['40_24'] * self._rate_res[39])

    def _dCdt(self, t, mo_comps, vol, flow, in_comps, fix_DO, DO_sat_T):
        '''
        Defines dC/dt for the reactor based on mass balance.
        
        Overall mass balance:
        dComp/dt == InfFlow / Actvol * (in_comps - mo_comps) + GrowthRate
                 == (in_comps - mo_comps) / HRT + GrowthRate
                 
        Args:
            t (float): Time (days).
            mo_comps (list): Mass of each component (mg/L).
            vol (float): Reactor volume (m3).
            flow (float): Influent flow rate (m3/d).
            in_comps (list): Influent component concentrations (mg/L).
            fix_DO (bool): Whether to fix the DO concentration.
            DO_sat_T (float): Saturation DO at the chosen temperature (mg/L).

        Returns:
            dCdt (list): List of dC/dt for each component (mg/L/d).

        '''

        # Load all normalized reaction rates
        self._reaction_rate(mo_comps)

        # Calculate hydraulic retention time
        _HRT = vol / flow

        if fix_DO or self._bulk_DO == 0:
            result = [0.0]
        else:
            result = [(in_comps[0] - mo_comps[0]) / _HRT
                      + self._KLa * (DO_sat_T - mo_comps[0])
                      + self._rate0_S_O2()]
            
        result.append((in_comps[1] - mo_comps[1]) / _HRT
                      + self._rate1_S_F())
        
        result.append((in_comps[2] - mo_comps[2]) / _HRT
                        + self._rate2_S_A())
        
        result.append((in_comps[3] - mo_comps[3]) / _HRT
                        + self._rate3_S_NH4())
        
        result.append((in_comps[4] - mo_comps[4]) / _HRT
                        + self._rate4_S_NH2OH())
        
        result.append((in_comps[5] - mo_comps[5]) / _HRT
                        + self._rate5_S_N2O())
        
        result.append((in_comps[6] - mo_comps[6]) / _HRT
                        + self._rate6_S_NO())
        
        result.append((in_comps[7] - mo_comps[7]) / _HRT
                        + self._rate7_S_NO2())
        
        result.append((in_comps[8] - mo_comps[8]) / _HRT
                        + self._rate8_S_NO3())
        
        result.append((in_comps[9] - mo_comps[9]) / _HRT
                        + self._rate9_S_PO4())
        
        result.append((in_comps[10] - mo_comps[10]) / _HRT
                        + self._rate10_S_I())
        
        result.append((in_comps[11] - mo_comps[11]) / _HRT
                        + self._rate11_S_ALK())
        
        result.append((in_comps[12] - mo_comps[12]) / _HRT
                        + self._rate12_S_N2())
        
        result.append((in_comps[13] - mo_comps[13]) / _HRT
                        + self._rate13_X_I())
        
        result.append((in_comps[14] - mo_comps[14]) / _HRT
                        + self._rate14_X_S())
        
        result.append((in_comps[15] - mo_comps[15]) / _HRT
                        + self._rate15_X_H())
        
        result.append((in_comps[16] - mo_comps[16]) / _HRT
                        + self._rate16_X_PAO())
        
        result.append((in_comps[17] - mo_comps[17]) / _HRT
                        + self._rate17_X_PP())
        
        result.append((in_comps[18] - mo_comps[18]) / _HRT
                        + self._rate18_X_PHA())
        
        result.append((in_comps[19] - mo_comps[19]) / _HRT
                        + self._rate19_X_AOB())
        
        result.append((in_comps[20] - mo_comps[20]) / _HRT
                        + self._rate20_X_NOB())
        
        result.append((in_comps[21] - mo_comps[21]) / _HRT
                        + self._rate21_X_TSS())
        
        result.append((in_comps[22] - mo_comps[22]) / _HRT
                        + self._rate22_X_MeOH())
        
        result.append((in_comps[23] - mo_comps[23]) / _HRT
                        + self._rate23_X_MeP())
        
        return result
            