from math import fabs, log, log10, exp
from exxonmobil_irms.workflows.fluid_property.rules import FluidProperty
from exxonmobil_irms.units.decoration import units, UnitSystemType, UnitType
from exxonmobil_irms.core.utils import Logger_error, Logger_info, Logger_warning

print_info = False

class FluidPropertyHGOR(FluidProperty):
    def __init__(self, case_or_project, Pmin=None, Pmax=None, min_visc=None,
                 npts=200, save=True, unit_system=UnitSystemType.storage,
                 clear_config=False):
        super(FluidPropertyHGOR,self).__init__(case_or_project, Pmin, Pmax, min_visc,
                                               npts, save, unit_system, clear_config)
        Logger_warning("Use the modified black oil to model high GOR fluid behavior")

    def computeSolutionGasOilRatio(self,regionNum):
        """
        Parameters
        ----------
        regionNum : Integer
        """
        # Vazquez and Beggs, 1980 (Default in EMPower)
        # Source, Protherm
        # Gas spec. gravity: 0.511 - 1.351
        # Oil gravity:      15.3 - 59.5 API
        # Pressure:         141.0 - 9515.0 psia
        # Soln GOR:         90.3 - 2199.0 Scf/Stb
        # Temperature:      100.0 - 258.0 F
        if self.FluidType[regionNum] == 'Deadoil':
            self.Rso[regionNum]  = [0.0 for _ in range(self.npts)]
        else:
            C1 = 0.0362
            C2 = 1.0937
            C3 = 25.7240
            if(self.API[regionNum] > (30.0 + 1e-12)):
                C1 = 0.0178
                C2 = 1.1870
                C3 = 23.9310
            Psat     = self.Psat[regionNum]
            Tsat     = self.Tres[regionNum] # assumed as Tsat
            API      = self.API[regionNum]
            Gamma_gs = self.Gamma_gs[regionNum]
            a        = C1*Gamma_gs
            c        = exp(C3 * API / (Tsat + 459.67))
            self.Rso[regionNum] = [a*(press**C2)*c for press in self.Pressure if press <= Psat]    #[(C1*Gamma_gs)*(press**C2)*exp(C3*API/(Tsat+459.67)) for press in self.Pressure if press <= Psat]
            self.extrapolateSolutionGORAbovePsat(regionNum)

    def extrapolateSolutionGORAbovePsat(self,regionNum):
        """
        Parameters
        ----------
        regionNum : Integer
        """
        # Simple linear extrapolation
        indx      = self.index_Psat[regionNum] - 1
        slope     = (self.Rso[regionNum][indx] - self.Rso[regionNum][indx-1])/(self.Pressure[indx] - self.Pressure[indx-1])
        if(self.UnderSaturatedOilExtrapolationMethod == 'Constant' or self.UnderSaturatedOilExtrapolationMethod == 'Default'):
            slope = 0.0
            if print_info:
                Logger_info("Constant GOR Extrapolation above Bubble Point")
        else:
            if print_info:
                Logger_info("Linear GOR Extrapolation above Bubble Point")

        intercept = self.Rso[regionNum][indx] - slope*self.Pressure[indx]
        for i in range(indx+1,self.npts):
            self.Rso[regionNum].append(slope*self.Pressure[i] + intercept)

    def computeLiveOilFVF(self,regionNum):
        """
        Parameters
        ----------
        regionNum : Integer
        """
        #Vasques and Beggs, 1980
        # Property Range
        # Gas spec. gravity: 0.511 - 1.351
        # Oil FVF: 1.02 - 3.0 Bbl/Stb
        # Oil gravity: 15.3 - 59.5 API
        # Soln GOR: 90.3 - 2199.0 Scf/Stb
        # Temperature: 100.0 - 258.0 F
        C1 = 4.677e-4
        C2 = 1.751e-5
        C3 =-1.811e-8
        if(self.API[regionNum] > (30.0 + 1e-12)):
            C1 = 4.670e-4
            C2 = 1.100e-5
            C3 = 1.337e-9

        Psat     = self.Psat[regionNum]
        Tres     = self.Tres[regionNum] # assumed as Tsat
        API      = self.API[regionNum]
        Gamma_gs = self.Gamma_gs[regionNum]
        C1       = C1*1.0
        C2       = C2*(Tres-60)*(API/Gamma_gs)
        C3       = C3*(Tres-60)*(API/Gamma_gs)
        self.Bo[regionNum] = [1.0 + C1*self.Rso[regionNum][i] + C2 + C3*self.Rso[regionNum][i] for i, press in enumerate(self.Pressure) if press <= Psat]

        self.extrapolateLiveOilFVFAbovePsat(regionNum)
        self.computeLiveOilDensity(regionNum)

    def computeIsothermalLiveOilCompressibilityAbovePsat(self, regionNum):
        """
        Parameters
        ----------
        regionNum : Integer
        """
        # Vasques and Beggs, 1980, From EMPower
        Tres     = self.Tres[regionNum]         # F
        API      = self.API[regionNum]
        Gamma_gs = self.Gamma_gs[regionNum]
        for indx in range(self.index_Psat[regionNum],self.npts):
            GOR  = self.Rso[regionNum][indx]
            self.Co[regionNum][indx] = (-1433.0 + 5.0*GOR + 17.2*Tres - 1180.0*Gamma_gs + 12.61*API)/(self.Pressure[indx]*1e+5)

        self.computeIsothermalLiveOilCompressibilityBelowPsat(regionNum)

    def computeIsothermalLiveOilCompressibilityBelowPsat(self, regionNum):
        """
        Parameters
        ----------
        regionNum : Integer
        """
        # Based on McCain and Coathurs, 1988 (refer Tareq Ahmed, Page 226 PVT Book
        for indx, press in enumerate(self.Pressure):
            if(press <= self.Psat[regionNum]):
                A = -7.573 - 1.45 * log(press) - 0.383 * log(self.Psat[regionNum]) + 1.402 * log(self.Tres[regionNum]) + 0.256 * log(self.API[regionNum]) + 0.449 * log(self.GOR[regionNum])
                self.Co[regionNum][indx] = exp(A)

    def extrapolateLiveOilFVFAbovePsat(self, regionNum):
        """
        Parameters
        ----------
        regionNum : Integer
        """
        # EMPower, 2000
        self.computeIsothermalLiveOilCompressibilityAbovePsat(regionNum)
        #indx      = self.index_Psat[regionNum]
        #Psat      = self.Psat[regionNum]
        #Bo        = self.Bo_sat[regionNum]
        #for i in range(indx,self.npts):
        #    Cbo = (1/exp(self.Co[regionNum][i] *(Psat - self.Pressure[i])) -1)/(self.Pressure[i]-Psat)
        #    self.Bo[regionNum].append(Bo/(1.0 + Cbo*(self.Pressure[i] - Psat)))
        #    #self.Bo[regionNum].append(Bo/(1.0 + self.Co[regionNum][i]*(self.Pressure[i] - Psat)))

        # Simple linear extrapolation
        indx      = self.index_Psat[regionNum] - 1
        slope     = (self.Bo[regionNum][indx] - self.Bo[regionNum][indx-1])/(self.Pressure[indx] - self.Pressure[indx-1])
        slope = slope/50.
        intercept = self.Bo[regionNum][indx] - slope*self.Pressure[indx]
        if(self.UnderSaturatedOilExtrapolationMethod == 'Linear'):
            for i in range(indx+1,self.npts):
                # simple Linear Extrapolation
                self.Bo[regionNum].append(slope*self.Pressure[i] + intercept)
            if print_info:
                Logger_info("Bo Extrapolation based on Simple Linear Extrapolation")
        elif self.UnderSaturatedOilExtrapolationMethod in ('Default', None):
            for i in range(indx+1,self.npts):
                # For constant compressiblity:
                self.Bo[regionNum].append(self.Bo_sat[regionNum] * exp(-self.Co[regionNum][i] * (self.Pressure[i] - self.Psat[regionNum])))
                #For compressibility evaluated based on Vasquez and Beggs, use the following:
                #A = 1e-5*(-1433 + 5*self.GOR[regionNum] + 17.2*self.Tres[regionNum] - 1180*Gamma_gs + 12.61*API)
                #self.Bo[regionNum].append(self.Bo_sat[regionNum]*exp(-A*log(self.Pressure[i]/self.Psat[regionNum])))
            if print_info:
                Logger_info("Bo Extrapolation based on Constant Compressibility Defintion")


    def computeLiveOilDensity(self, regionNum):
        """
        Parameters
        ----------
        regionNum : Integer
        """
        oilSp    = self.OilSpGravity[regionNum]
        gasSp    = self.GasSpGravity[regionNum]

        for indx, press in enumerate(self.Pressure):
            if press <= self.Psat[regionNum]:
                #from EMPower correlation
                #self.OilDensity[regionNum][indx] = (350*oilSp + 0.0764*gasSp*self.Rso[regionNum][indx])/(5.6146*self.Bo[regionNum][indx]*62.428) # in gm/cc

                # From material balance: (Note check multiplier in front of Oil Sp)
                self.OilDensity[regionNum][indx]= (62.428*oilSp + 0.0136*self.Rso[regionNum][indx]*gasSp)/(62.428*self.Bo[regionNum][indx]) # in gm/cc

                # if using vazquez and beggs correlation for compressibility, use the following eqt for density
                #A = 1e-5*(-1433 + 5*self.Rso[regionNum][indx] + 17.2*self.Tres[regionNum] - 1180*Gamma_gs* + 12.61*API)
                #self.OilDensity[regionNum][indx] =  self.OilDensity_sat[regionNum]*exp(A*log((press/self.Psat[regionNum])))/62.428 # in g/cc

        # below equation is EMPower's formulation -> for constant compressibility
        for indx in range(self.index_Psat[regionNum],self.npts):
            press = self.Pressure[indx]
            self.OilDensity[regionNum][indx] = (self.OilDensity_sat[regionNum] * exp(self.Co[regionNum][indx] * (press - self.Psat[regionNum]))) # in gm/cc

            # Material Balance
            #self.OilDensity[regionNum][indx] = (62.4*oilSp + 0.0136*self.Rso[regionNum][indx]*gasSp)/self.Bo[regionNum][indx]

    def computeLiveOilViscosity(self,regionNum):
        """
        Parameters
        ----------
        regionNum : Integer
        """
        # Beggs and Robinson, 1975 Defailt EMPower
        self.computeDeadOilViscosity(regionNum)
        for indx in range(0,self.index_Psat[regionNum]):
            a             = 10.715*((self.Rso[regionNum][indx] + 100.0)**(-0.515))
            b             = 5.44*  ((self.Rso[regionNum][indx] + 150.0)**(-0.338))
            self.Visc_oil[regionNum][indx] = a*(self.Visc_oil[regionNum][indx]**b)
        self.computeLiveOilViscosityAboveBubblePt(regionNum)

    def computeDeadOilViscosity(self, regionNum):
        """
        Parameters
        ----------
        regionNum : Integer
        """
        # Beggs and Robinson, 1975 Default EMPower
        b             = exp(6.9824 - 0.04658 * self.API[regionNum])
        a             = b*(self.Tres[regionNum]**(-1.163))
        self.Visc_oil[regionNum] = [10**a - 1]*self.npts

    def computeLiveOilViscosityAboveBubblePt(self, regionNum):
        """
        Parameters
        ----------
        regionNum : Integer
        """
        if(self.UnderSaturatedOilExtrapolationMethod == 'Default'):
            # Vazques and Beggs, 1980
            for indx in range(self.index_Psat[regionNum],self.npts):
                # Below is EMPower formulation
                m    = 2.6 * (self.Pressure[indx]**1.187) * exp(-11.513 - 8.98e-5 * self.Pressure[indx])
                # Fomrulation given in Tariq
                #A     =-3.9*(1e-5)*self.Pressure[indx] - 5
                #m     = 2.6*(self.Pressure[indx]**1.187)*(10**A)
                self.Visc_oil[regionNum][indx] = self.Visc_oil_sat[regionNum]*(self.Pressure[indx]*1.0/self.Psat[regionNum])**m
        else:
            indx      = self.index_Psat[regionNum] - 1
            slope     = (self.Visc_oil[regionNum][indx] - self.Visc_oil[regionNum][indx-1])/(self.Pressure[indx] - self.Pressure[indx-1])
            intercept = self.Visc_oil[regionNum][indx] - slope*self.Pressure[indx]
            if (slope*self.Pmax + intercept ) < 0.0:
                min_viscosity = self.min_viscosity
                slope = (min_viscosity - self.Visc_oil[regionNum][indx])/(self.Pmax - self.Pressure[indx])
                intercept = self.Visc_oil[regionNum][indx] - slope*self.Pressure[indx]
            for indx in range(self.index_Psat[regionNum],self.npts):
                self.Visc_oil[regionNum][indx] = slope*self.Pressure[indx] + intercept
            if print_info:
                Logger_info("Linear Viscosity Extrapolation above Bubble Point")

    def computeDryGasFVF(self, regionNum):
        """
        Parameters
        ----------
        regionNum : Integer
        """
        self.computeDryGasZFactor(regionNum)
        fac = self.Pstd[regionNum]/(self.Tstd[regionNum] + 459.67) # use this fro ft3/scf
        fac = fac*0.178107607 # conversion factor for cubic feet to bbl for us crude oil http://www.asknumbers.com/CubicFeetToBarrel.aspx
        for indx, press in enumerate(self.Pressure):
            self.Bg[regionNum][indx] = fac*self.Zfactor[regionNum][indx]*(self.Tres[regionNum] + 459.67)/press

    def computeDryGasZFactor(self,regionNum):
        """
        Parameters
        ----------
        regionNum : Integer
        """
        # Dranchuk and Abou-Kassem, 1975-Default EMPower
        Zr                      = 1.0
        self.Zfactor[regionNum] = [Zr]*self.npts
        A1       = 0.3265
        A2       =-1.0700
        A3       =-0.5339
        A4       = 0.01569
        A5       =-0.05165
        A6       = 0.5475
        A7       =-0.7361
        A8       = 0.1844
        A9       = 0.1056
        A10      = 0.6134
        A11      = 0.7210
        for indx, press in enumerate(self.Pressure):
            error= self.LARGE
            iter = 0
            Tpr  = (self.Tres[regionNum] + 459.67)/self.Tpc[regionNum]
            Ppr  = press/self.Ppc[regionNum]
            Rpr  = 0.27*Ppr/(self.Zfactor[regionNum][indx]*Tpr)
            assert Tpr >= 1.0 and Tpr <= 3.0,  'Pseudo Reduced Temperature, Tpr: ' + str(Tpr) + ' for Region: ' + str(regionNum) + ' is Out Of Bounds: 1.0 <= Tpr <=3.0'
            #assert Ppr >= 0.2 and Ppr <= 30.0, 'Pseudo Reduced Pressure   , Ppr: ' + str(Ppr) + ' for Region: ' + str(regionNum) + ' is Out Of Bounds: 0.2 <= Ppr <=30.0'
            # Newton Raphson to evaluate Density
            while(error > self.TINY and iter < self.iterMax):
                #Note: Tpc is computed in Rankine according to correlation
                #Starling-Carnahan equation of state
                # 1.0 <= Tpr <=3.0
                # 0.2 <= Ppr <= 30.0 and
                Rpr_Old = Rpr
                a       = (A1 + A2/Tpr + A3/(Tpr**3) + A4/(Tpr**4) + A5/(Tpr**5))
                b       = 0.27*Ppr/Tpr
                c       = A6 + A7/Tpr + A8/(Tpr**2)
                d       = A9*(A7/Tpr + A8/(Tpr**2))
                e       = A10/(Tpr**3)
                Zr      = 1.0 + a * Rpr - b / Rpr + c * (Rpr**2) - d * (Rpr**5) + e * (1 + A11*(Rpr**2)) * (Rpr**2) * exp(-A11 * (Rpr ** 2))
                Zprime  = a + b / (Rpr**2) + 2 * c * Rpr - 5 * d * (Rpr**4) + 2 * e * Rpr * exp(-A11 * (Rpr ** 2)) * (1 + 2 * A11 * (Rpr ** 3) - A11 * (Rpr ** 2) * (1 + A11 * (Rpr ** 2)))
                Rpr     = Rpr_Old - Zr/Zprime
                error   = fabs(Rpr - Rpr_Old)
                iter   += 1
            Zr = 0.27*Ppr/(Rpr*Tpr)
            assert Zr > 0.0, 'Error in Compressibility Computation at indx: ' + str(indx) + ' region ID: ' + str(regionNum) + ' error: ' + str(error) + ' iter: ' + str(iter) + ' Zr: ' + str(Zr)
            self.Zfactor[regionNum][indx] = Zr

        self.computeGasCompressibility(regionNum)

    def computeGasCompressibility(self, regionNum):
        """
        Parameters
        ----------
        regionNum : Integer
        """
        A1       = 0.3265
        A2       =-1.0700
        A3       =-0.5339
        A4       = 0.01569
        A5       =-0.05165
        A6       = 0.5475
        A7       =-0.7361
        A8       = 0.1844
        A9       = 0.1056
        A10      = 0.6134
        A11      = 0.7210
        for indx, press in enumerate(self.Pressure):
            Tpr     = (self.Tres[regionNum] + 459.67)/self.Tpc[regionNum]
            Ppr     = press/self.Ppc[regionNum]
            Rpr     = 0.27*Ppr/(self.Zfactor[regionNum][indx]*Tpr)
            a       = (A1 + A2/Tpr + A3/(Tpr**3) + A4/(Tpr**4) + A5/(Tpr**5))
            c       = A6 + A7/Tpr + A8/(Tpr**2)
            d       = A9*(A7/Tpr + A8/(Tpr**2))
            e       = A10/(Tpr**3)
            Zr      = self.Zfactor[regionNum][indx]
            Zprime  = a + 2 * c * Rpr - 5 * d * (Rpr**4) + 2 * e * Rpr * exp(-A11 * (Rpr ** 2)) * (1 + 2 * A11 * (Rpr ** 3) - A11 * (Rpr ** 2) * (1 + A11 * (Rpr ** 2)))
            dZdPpr  = (0.27/(Zr*Zr*Tpr))*(Zprime/(1.0 + (Rpr/Zr)*Zprime))
            Cpr     = 1.0/Ppr - (1.0/Zr)*dZdPpr
            assert Cpr > 0.0, 'Error in Compressibility Computation at indx: ' + str(indx) + ' region ID: ' + str(regionNum) + ' Zr: ' + str(Zr) + ' dZdPpr: ' + str(dZdPpr) + ' Cpr: ' + str(Cpr)
            self.Cg[regionNum][indx] = Cpr/self.Ppc[regionNum]

    def computeDryGasViscosity(self, regionNum):
        """
        Parameters
        ----------
        regionNum : Integer
        """
        self.computeDryGasDensity(regionNum)
        Mg = self.GasMa[regionNum]
        A = (9.379 + 0.01607*Mg)*(self.Tres[regionNum]+ 459.67)**1.5/(209.2 + 19.26*Mg + (self.Tres[regionNum]+ 459.67))
        B = 3.448 + (986.4/(self.Tres[regionNum]+ 459.67)) + 0.01009*Mg
        C = 2.447 - 0.2224*B
        self.Visc_gas[regionNum] = [A * 1e-4 * exp(B * (den ** C)) for den in self.GasDensity[regionNum]]

    def computeDryGasDensity(self, regionNum):
        """
        Parameters
        ----------
        regionNum : Integer
        """
        fac = self.GasMa[regionNum]/self.GasConstant*self.lbFt3ToGmCc
        # Den = Ma*P/ZRT
        self.GasDensity[regionNum] = [fac*press/(self.Zfactor[regionNum][indx]*(self.Tres[regionNum] + self.FarToRankine)) for indx, press in enumerate(self.Pressure)]


    def computerGasWaterRatio(self, regionNum):
        """
        Parameters
        ----------
        regionNum : Integer
        """
        Tr= self.Tres[regionNum]
        # Weight percent dissolved solids, Cw (Cw = Cppm *10-4)
        S     = self.Salinity[regionNum]*1e-4
        A     =  8.15839      - 6.12265*1e-2*Tr + 1.91663*1e-4*(Tr**2) - 2.1654 *1e-7*(Tr**3)
        B     =  1.01021*1e-2 - 7.44241*1e-5*Tr + 3.05553*1e-7*(Tr**2) - 2.94883*1e-10*(Tr**3)
        C     = -1e-7*(9.02505- 0.130237*Tr     + 8.53425*1e-4*(Tr**2) - 2.34122*1e-6*(Tr**3)+ 2.37049*1e-9*(Tr**4))
        D     = -0.0840655*S*(Tr**(-0.285854))
        Pctrl = self.Pres[regionNum]
        if self.Psat[regionNum] is not None:
            Pctrl = self.Psat[regionNum]

        if self.FluidType[regionNum] == 'Drygas':
            Pctrl = self.LARGE

        for indx, press in enumerate(self.Pressure):
            if press <= Pctrl:
                Rsw                       = A + B*press + C*(press**2)
                self.Rsw[regionNum][indx] = Rsw*(10**D)

        if self.FluidType[regionNum] != 'Drygas':
            self.extrapolateGasWaterRatioAboveBubblePt(regionNum)

    def extrapolateGasWaterRatioAboveBubblePt(self, regionNum):
        """
        Parameters
        ----------
        regionNum : Integer
        """
        indx      = self.index_Psat[regionNum]-1
        slope     =(self.Rsw[regionNum][indx] - self.Rsw[regionNum][indx-1])/(self.Pressure[indx] - self.Pressure[indx-1])
        intercept = self.Rsw[regionNum][indx] - slope*self.Pressure[indx]
        for i in range(indx,self.npts):
            self.Rsw[regionNum][i] = (slope*self.Pressure[i] + intercept)

    def computeWaterFVF(self, regionNum):
        """
        Parameters
        ----------
        regionNum : Integer
        """
        dVwt  = -1.0001e-2 + 1.33391e-4*self.Tres[regionNum] + 5.50654e-7*(self.Tres[regionNum]**2)
        Pctrl = self.Pres[regionNum]

        if self.Psat[regionNum] is not None:
            Pctrl = self.Psat[regionNum]

        if self.FluidType[regionNum] == 'Drygas':
            Pctrl = self.LARGE

        for indx, press in enumerate(self.Pressure):
            if press <= Pctrl:
                dVwp = -1.95301e-9*self.Tres[regionNum]*press - 1.72834e-13*self.Tres[regionNum]*(press**2) - 3.58922e-7*press - 2.25341e-10*(press**2)
                self.Bw[regionNum][indx] = (1+dVwt)*(1+dVwp)

        self.computeIsothermalWaterCompressiblity(regionNum)

        if self.FluidType[regionNum] != 'Drygas':
            self.computeWaterFVFAboveBubblePt(regionNum)

        self.computeWaterDensity(regionNum)

    def computeWaterFVFAboveBubblePt(self, regionNum):
        """
        Parameters
        ----------
        regionNum : Integer
        """
        press = self.Pres[regionNum] # reservoir pressure and not self.Pressure
        if self.Psat[regionNum] is not None:
            press = self.Psat[regionNum]
        # Compute water FVF at bubble point
        dVwt  = -1.0001e-2 + 1.33391e-4*self.Tres[regionNum] + 5.50654e-7*(self.Tres[regionNum]**2)
        dVwp  = -1.95301e-9*self.Tres[regionNum]*press - 1.72834e-13*self.Tres[regionNum]*(press**2) - 3.58922e-7*press - 2.25341e-10*(press**2)
        Bwb   = (1.0+dVwt)*(1.0+dVwp)
        for indx, pr in enumerate(self.Pressure):
            if pr > press:
                self.Bw[regionNum][indx] = Bwb /(1.0 + self.Cw[regionNum][indx]*(pr -press))

    def computeWaterDensity(self, regionNum):
        """
        Parameters
        ----------
        regionNum : Integer
        """
        for indx, Bw in enumerate(self.Bw[regionNum]):
            self.WaterDensity[regionNum][indx] = self.WaterDensity_sc[regionNum]/(Bw*(1.0 + (self.Rsw[regionNum][indx]*self.GasDensity_sc[regionNum])/(5.6146*self.WaterDensity_sc[regionNum])))

    def computeIsothermalWaterCompressiblity(self, regionNum):
        """
        Parameters
        ----------
        regionNum : Integer
        """
         #cw = Isothermal water (brine) compressibility, psi-1
         #Pr = Reservoir pressure, psia
         #S = Salinity, mg/L, Cs (Cs = denw(sc) * Cppm), Cppm = Dissolved solids parts per million (Cppm = Cw *10^4)
         #Tr = Reservoir temperature, F
         #denw(sc) = Water density at standard conditions, gm/cc
        # Osif, 1988 - Default
        salinity = self.Salinity[regionNum]*self.WaterDensity_sc[regionNum]
        for indx, press in enumerate(self.Pressure):
            self.Cw[regionNum][indx] = 1.0/(7.033*press + 0.5415*salinity - 537.0*self.Tres[regionNum] + 403300.0)

    def computerWaterViscosity(self,regionNum):
        """
        Parameters
        ----------
        regionNum : Integer
        """
        for indx, press in enumerate(self.Pressure):
            self.Visc_water[regionNum][indx] = self.Visc_water_sc[regionNum]*(0.9994 + 4.0295*1e-5*press + 3.1062*1e-9*(press**2))

