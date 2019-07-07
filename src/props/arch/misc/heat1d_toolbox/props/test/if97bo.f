C     Last change:  RW   24 Mar 2003    2:54 pm
C*************************************************************************
C                                                                        *
C                       RUHR-UNIVERSITAET BOCHUM                         *
C                     Fakultaet fuer Maschinenbau                        *
C                     Lehrstuhl fuer Thermodynamik                       *
C                       Prof. Dr.-Ing. W. Wagner                         *
C                           D - 44780 Bochum                             *
C                                                                        *
C                         Tel.: +49-234-32 23033                         *
C                         Fax:  +49-234-32 14163                         *
C                                                                        *
C                                                                        *
C*************************************************************************
C
C  SOFTWARE FOR THE IMPLEMENTATION OF THE
C
C  "IAPWS INDUSTRIAL FORMULATION FOR THE THERMODYNAMIC PROPERTIES
C  OF WATER AND STEAM (IAPWS-IF97)"
C
C  A DESCRIPTION OF THE FUNCTIONS IS GIVEN IN THE FILE README.TXT.
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION VBPT(P,T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF SPECIFIC VOLUME
C FOR GIVEN VALUES OF P AND T
C
C INPUT:   P  Pressure [MPa]
C          T  Temperature [Kelvin]
C
C OUTPUT:  V     Specific Volume [cubic meter per kilogram]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSOPT(P,T,IREG)
C  
      IF (IREG .NE. 0) THEN
         VBPT=VPT(P,T,IREG)
      ELSE
         VBPT=0.D0
      ENDIF   
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION DBPT(P,T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF DENSITY
C FOR GIVEN VALUES OF P AND T
C
C INPUT:   P  Pressure [MPa]
C          T  Temperature [Kelvin]
C
C OUTPUT:  D     Density [kilogram per cubic meter]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSOPT(P,T,IREG)
C  
      IF (IREG .NE. 0) THEN
         DBPT=1.D0/VPT(P,T,IREG)
      ELSE
         DBPT=0.D0
      ENDIF   
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION HBPT(P,T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF SPECIFIC ENTHALPY
C FOR GIVEN VALUES OF P AND T
C
C INPUT:   P  Pressure [MPa]
C          T  Temperature [Kelvin]
C
C OUTPUT:  H     Specific Enthalpy [kilojoule per kilogram]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSOPT(P,T,IREG)
C  
      IF (IREG .NE. 0) THEN
         HBPT=HPT(P,T,IREG)
      ELSE   
         HBPT=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION SBPT(P,T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF SPECIFIC ENTROPY
C FOR GIVEN VALUES OF P AND T
C
C INPUT:   P  Pressure [MPa]
C          T  Temperature [Kelvin]
C
C OUTPUT:  S     Specific Entropy [kilojoule per kilogram and Kelvin]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSOPT(P,T,IREG)
C  
      IF (IREG .NE. 0) THEN
         SBPT=SPT(P,T,IREG)
      ELSE   
         SBPT=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION CPBPT(P,T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF SPECIFIC ISOBARIC HEAT CAPACITY
C FOR GIVEN VALUES OF P AND T
C
C INPUT:   P  Pressure [MPa]
C          T  Temperature [Kelvin]
C
C OUTPUT:  CP    Specific Isobaric Heat Capacity 
C                [kilojoule per kilogram and Kelvin]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSOPT(P,T,IREG)
C  
      IF (IREG .NE. 0) THEN
         CPBPT=CPPT(P,T,IREG)
      ELSE
         CPBPT=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION CVBPT(P,T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF SPECIFIC ISOCHORIC HEAT CAPACITY
C FOR GIVEN VALUES OF P AND T
C
C INPUT:   P  Pressure [MPa]
C          T  Temperature [Kelvin]
C
C OUTPUT:  CV    Specific Isochoric Heat Capacity 
C                [kilojoule per kilogram and Kelvin]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSOPT(P,T,IREG)
C  
      IF (IREG .NE. 0) THEN
         CVBPT=CVPT(P,T,IREG)
      ELSE
         CVBPT=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION WBPT(P,T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF SPEED OF SOUND
C FOR GIVEN VALUES OF P AND T
C
C INPUT:   P  Pressure [MPa]
C          T  Temperature [Kelvin]
C
C OUTPUT:  W     Speed of sound [meter per second]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSOPT(P,T,IREG)
C  
      IF (IREG .NE. 0) THEN
         WBPT=WPT(P,T,IREG)
      ELSE
         WBPT=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION UBPT(P,T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF SPECIFIC INTERNAL ENERGY
C FOR GIVEN VALUES OF P AND T
C
C INPUT:   P  Pressure [MPa]
C          T  Temperature [Kelvin]
C
C OUTPUT:  U     Specific internal energy [kilojoule per kilogram]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSOPT(P,T,IREG)
C  
      IF ((IREG .EQ. 1) .OR. (IREG .EQ. 2) .OR.
     & (IREG .EQ. 3) .OR. (IREG .EQ. 5)) THEN
         UBPT=UPT(P,T,IREG)
      ELSE
         UBPT=0.D0
      ENDIF   
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION FBPT(P,T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF SPECIFIC HELMHOLTZ FREE ENERGY
C FOR GIVEN VALUES OF P AND T
C
C INPUT:   P  Pressure [MPa]
C          T  Temperature [Kelvin]
C
C OUTPUT:  F     Specific Helmholtz free energy [kilojoule per kilogram]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSOPT(P,T,IREG)
C  
      IF ((IREG .EQ. 1) .OR. (IREG .EQ. 2) .OR.
     & (IREG .EQ. 3) .OR. (IREG .EQ. 5)) THEN
         U=UPT(P,T,IREG)
         S=SPT(P,T,IREG)
         FBPT=U-T*S
      ELSE
         FBPT=0.D0
      ENDIF   
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION GBPT(P,T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF SPECIFIC GIBBS FREE ENERGY
C FOR GIVEN VALUES OF P AND T
C
C INPUT:   P  Pressure [MPa]
C          T  Temperature [Kelvin]
C
C OUTPUT:  G     Specific Gibbs free energy [kilojoule per kilogram]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSOPT(P,T,IREG)
C  
      IF ((IREG .EQ. 1) .OR. (IREG .EQ. 2) .OR.
     & (IREG .EQ. 3) .OR. (IREG .EQ. 5)) THEN
         GBPT=GPT(P,T,IREG)
      ELSE
         GBPT=0.D0
      ENDIF   
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION CAPBPT(P,T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF ISENTROPIC EXPONENT
C FOR GIVEN VALUES OF P AND T
C
C INPUT:   P  Pressure [MPa]
C          T  Temperature [Kelvin]
C
C OUTPUT:  CAP   Isentropic Exponent [-]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSOPT(P,T,IREG)
C  
      IF ((IREG .EQ. 1) .OR. (IREG .EQ. 2) .OR.
     & (IREG .EQ. 3) .OR. (IREG .EQ. 5)) THEN
         V=VPT(P,T,IREG)
         W=WPT(P,T,IREG)
         CAPBPT=W*W/P/V*1.D-6
      ELSE
         CAPBPT=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION ALPBPT(P,T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF ISOBARIC VOLUME EXPANSION COEFFICIENT
C FOR GIVEN VALUES OF P AND T
C
C INPUT:   P  Pressure [MPa]
C          T  Temperature [Kelvin]
C
C OUTPUT:  ALP   Isobaric Volume Expansion Coefficient [1 / Kelvin]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSOPT(P,T,IREG)
C  
      IF (IREG .NE. 0) THEN
         V=VPT(P,T,IREG)
         ALPBPT=DVDTPT(P,V,T,IREG)/V
      ELSE
         ALPBPT=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION BETBPT(P,T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF ISOCHORIC PRESSURE COEFFICIENT
C FOR GIVEN VALUES OF P AND T
C
C INPUT:   P  Pressure [MPa]
C          T  Temperature [Kelvin]
C
C OUTPUT:  BET   Isochoric Pressure Coefficient [1 / Kelvin]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSOPT(P,T,IREG)
C  
      IF (IREG .NE. 0) THEN
         V=VPT(P,T,IREG)
         BETBPT=(-1.D0)*DVDTPT(P,V,T,IREG)/DVDPPT(P,V,T,IREG)/P
      ELSE
         BETBPT=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION GAMBPT(P,T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF ISOTHERMAL COMPRESSIBILITY COEFFICIENT
C FOR GIVEN VALUES OF P AND T
C
C INPUT:   P  Pressure [MPa]
C          T  Temperature [Kelvin]
C
C OUTPUT:  GAM   Isothermal Compressibility Coefficient [1 / bar]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSOPT(P,T,IREG)
C  
      IF (IREG .NE. 0) THEN
         V=VPT(P,T,IREG)
         GAMBPT=DVDPPT(P,V,T,IREG)/V*(-1.D0)
      ELSE
         GAMBPT=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION FUGBPT(P,T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF FUGACITY
C FOR GIVEN VALUES OF P AND T
C
C INPUT:   P  Pressure [MPa]
C          T  Temperature [Kelvin]
C
C OUTPUT:  FUG   Fagacity [MPa]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      PARAMETER (r= 0.461526D0)
C      
      CALL REGSOPT(P,T,IREG)
C  
      IF ((IREG .EQ. 1) .OR. (IREG .EQ. 2) .OR.
     & (IREG .EQ. 3) .OR. (IREG .EQ. 5)) THEN
         G=GPT(P,T,IREG)
         G0=G0PT(P,T)
         FUGBPT=P*DEXP((G-G0)/R/T)
      ELSE
         FUGBPT=0.D0
      ENDIF
C
      END
C***********************************************************************
      DOUBLE PRECISION FUNCTION ETABPT(P,T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF DYNAMIC VISCOSITY
C FOR GIVEN VALUES OF P AND T
C
C INPUT:   P  Pressure [MPa]
C          T  Temperature [Kelvin]
C
C OUTPUT:  ETA   Dynamic Viscosity [kilogram per meter and second]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSOPT(P,T,IREG)
C  
      IF ((IREG .NE. 0) .AND. (T .LE. 1173.15D+0)) THEN
         V=VPT(P,T,IREG)
         CALL VISCVT(P,V,T,IREG,ETA)
         ETABPT=ETA*1.D6
      ELSE
         ETABPT=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION VNUBPT(P,T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF KINEMATIC VISCOSITY
C FOR GIVEN VALUES OF P AND T
C
C INPUT:   P  Pressure [MPa]
C          T  Temperature [Kelvin]
C
C OUTPUT:  NY    Kinematic Viscosity [square meter per second]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSOPT(P,T,IREG)
C  
      IF ((IREG .NE. 0) .AND. (T .LE. 1173.15D0)) THEN
         V=VPT(P,T,IREG)
         CALL VISCVT(P,V,T,IREG,ETA)
         VNUBPT=ETA*V*1.D6
      ELSE
         VNUBPT=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION PRBPT(P,T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF PRANDTL NUMBER
C FOR GIVEN VALUES OF P AND T
C
C INPUT:   P  Pressure [MPa]
C          T  Temperature [Kelvin]
C
C OUTPUT:  PR    Prandtl Number [-]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSOPT(P,T,IREG)
C  
      IF ((IREG .NE. 0) .AND. (IREG .NE. 5)) THEN
         V=VPT(P,T,IREG)
         CALL VISCVT(P,V,T,IREG,ETA)
         CALL THCONVT(P,V,T,IREG,TCO)
         CP=CPPT(P,T,IREG)
         PRBPT=ETA*CP/TCO*1.D3
      ELSE
         PRBPT=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION TCOBPT(P,T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF THERMAL CONDUCTIVITY
C FOR GIVEN VALUES OF P AND T
C
C INPUT:   P  Pressure [MPa]
C          T  Temperature [Kelvin]
C
C OUTPUT:  TCO   Thermal Conductivity [Watt per meter and Kelvin]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSOPT(P,T,IREG)
C  
      IF ((IREG .NE. 0) .AND. (IREG .NE. 5)) THEN
         V=VPT(P,T,IREG)
         CALL THCONVT(P,V,T,IREG,TCO)
         TCOBPT=TCO*1.D3
      ELSE
         TCOBPT=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION EPSBPT(P,T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF STATIC DIELECTRIC CONSTANT
C FOR GIVEN VALUES OF P AND T
C
C INPUT:   P  Pressure [MPa]
C          T  Temperature [Kelvin]
C
C OUTPUT:  EPS   Static Dielectric Constant [-]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSOPT(P,T,IREG)
C  
      IF ((IREG .NE. 0) .AND. (T .LE. 1073.15D0)) THEN
         D=1.D0/VPT(P,T,IREG)
         EPSBPT=DIELC(T,D)
      ELSE
         EPSBPT=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION REFBPT(P,T,WAV)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF REFRACTIVE INDEX
C FOR GIVEN VALUES OF P, T and WAV
C
C INPUT:   P  Pressure [MPa]
C          T  Temperature [Kelvin]
C          WAV    Wavelength of Light [meter * 1.D-6]
C
C OUTPUT:  REF   Refractive Index [-]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSOPT(P,T,IREG)
C  
      IF ((IREG .NE. 0) .AND. (T .LE. 773.15D0) .AND.
     & (WAV .GE. 0.2D0) .AND. (WAV .LE. 1.1D0)) THEN
         D=1.D0/VPT(P,T,IREG)
         REFBPT=REFI(T,D,WAV)
      ELSE
         REFBPT=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION DDDHBPT(P,T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF THE DERIVATIVE OF DENSITY WITH RESPECT 
C TO ENTHALPY AT CONSTANT PRESSURE
C FOR GIVEN VALUES OF P AND T
C
C INPUT:   P  Pressure [MPa]
C          T  Temperature [Kelvin]
C
C OUTPUT:  DDDH  Derivative of Density with Respect to Enthalpy at 
C                Constant Pressure 
C                [square kilogram per cubic meter and kilojoule]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSOPT(P,T,IREG)
C  
      IF (IREG .NE. 0) THEN
         V=VPT(P,T,IREG)
         CP=CPPT(P,T,IREG)
         DDDHBPT=(-1.D0)*DVDTPT(P,V,T,IREG)/CP/V/V
      ELSE
         DDDHBPT=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION DVDHBPT(P,T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF DERIVATIVE OF SPECIFIC VOLUME WITH 
C RESPECT TO ENTHALPY AT CONSTANT PRESSURE 
C FOR GIVEN VALUES OF P AND T
C
C INPUT:   P  Pressure [MPa]
C          T  Temperature [Kelvin]
C
C OUTPUT:  DVDH  Derivative of Specific Volume with respect to Enthalpy 
C                at Constant Pressure [cubic meter per kilojoule]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSOPT(P,T,IREG)
C  
      IF (IREG .NE. 0) THEN
         V=VPT(P,T,IREG)
         CP=CPPT(P,T,IREG)
         DVDHBPT=DVDTPT(P,V,T,IREG)/CP
      ELSE
         DVDHBPT=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION DDDPBPT(P,T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF DERIVATIVE OF DENSITY WITH RESPECT TO 
C PRESSURE AT CONSTANT ENTHALPY
C FOR GIVEN VALUES OF P AND T
C
C INPUT:   P  Pressure [MPa]
C          T  Temperature [Kelvin]
C
C OUTPUT:  DDDP  Derivative of Density with Respect to Pressure at 
C                Constant Enthalpy 
C                [kilogram per cubic meter and bar]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSOPT(P,T,IREG)
C  
      IF (IREG .NE. 0) THEN
         V=VPT(P,T,IREG)
         W=WPT(P,T,IREG)
         CP=CPPT(P,T,IREG)
         DDDPBPT=(DVDTPT(P,V,T,IREG)/CP/V+1.D0/W/W*1.D3)*1.D2
      ELSE
         DDDPBPT=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION DVDPBPT(P,T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF DERIVATIVE OF SPECIFIC VOLUME WITH 
C RESPECT TO PRESSURE AT CONSTANT ENTHALPY
C FOR GIVEN VALUES OF P AND T
C
C INPUT:   P  Pressure [MPa]
C          T  Temperature [Kelvin]
C
C OUTPUT:  DVDP  Derivative of Specific Volume with Respect to Pressure 
C                at Constant Enthalpy 
C                [cubic meter per kilogram and bar]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSOPT(P,T,IREG)
C  
      IF (IREG .NE. 0) THEN
         V=VPT(P,T,IREG)
         W=WPT(P,T,IREG)
         CP=CPPT(P,T,IREG)
         DVDPBPT=(DVDTPT(P,V,T,IREG)/CP+V/W/W*1.d3)*V*(-1.D2)
      ELSE
         DVDPBPT=0.D0
      ENDIF
C
      END
C
C*************************************************************************************
      DOUBLE PRECISION FUNCTION DHDPBPT(P,T)
C*************************************************************************************
C
C FUNCTION FOR THE CALCULATION OF DERIVATICVE OF ENTHALPY WITH RESPECT 
C TO PRESSURE AT CONSTANT TEMPERATURE
C FOR GIVEN VALUES OF P AND T
C
C INPUT:   P  Pressure [MPa]
C          T  Temperature [Kelvin]
C
C OUTPUT:  DHDP  Derivative of Enthalpy with Respect to Pressure at 
C                Constant Temperature 
C                [kilojoule per kilogram and bar]
C
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      CALL REGSOPT(P,T,IREG)
C
      IF (IREG .EQ. 1) THEN
           DHDPBPT=DHDP1(P,T)*1.D-1
      ELSEIF (IREG .EQ. 2) THEN
           DHDPBPT=DHDP2(P,T)*1.D-1
      ELSEIF (IREG .EQ. 3) THEN
           V=VPT3N(P,T)
           DHDV=DHDV3(V,T)
           DPDV=DPDV3(V,T)
           DHDPBPT=DHDV/DPDV*1.D-1
      ELSEIF (IREG .EQ. 5) THEN
           DHDPBPT=DHDP5(P,T)*1.D-1
      ELSE
           DHDPBPT=0.D0
      ENDIF
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION SIGBT(T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF SURFACE TENSION
C FOR GIVEN VALUES OF P AND T
C
C INPUT:   T  Temperature [Kelvin]
C
C OUTPUT:  SIG   Surface Tension [1.D-3 Newton per meter]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      IF ((T .GE. 273.15D0) .AND. (T .LE. 647.096D0)) THEN
         SIGBT=SURFTT(T)
      ELSE
         SIGBT=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION TBPH(P,H)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF TEMPERATURE
C FOR GIVEN VALUES OF P AND H
C
C INPUT:   P  Pressure [MPa]
C          H     Specific Enthalpy [kilojoule per kilogram]
C
C OUTPUT:  T     Temperature [Kelvin]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSPH(P,H,IREG)
C  
      IF ((IREG .EQ. 1) .OR. (IREG .EQ. 2) .OR. 
     1    (IREG .EQ. 3) .OR. (IREG .EQ. 5)) THEN
         TBPH=TPH(P,H,IREG)
      ELSEIF (IREG .EQ. 9) THEN
         TBPH=TSATPN(P)
      ELSE
         TBPH=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION VBPH(P,H)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF SPECIFIC VOLUME
C FOR GIVEN VALUES OF P AND H
C
C INPUT:   P  Pressure [MPa]
C          H     Specific Enthalpy [kilojoule per kilogram]
C
C OUTPUT:  V     Specific Volume [cubic meter per kilogram]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSPH(P,H,IREG)
C
      IF ((IREG .EQ. 1) .OR. (IREG .EQ. 2) .OR. 
     1    (IREG .EQ. 3) .OR. (IREG .EQ. 5)) THEN
         TB=TPH(P,H,IREG)
         VBPH=VPT(P,TB,IREG)
      ELSEIF (IREG .EQ. 9) THEN
         TS=TSATPN(P)
         IF (TS .LE. 623.15D0) THEN
            V1=VPT(P,TS,1)
            V2=VPT(P,TS,2)
            H1=HPT(P,TS,1)
            H2=HPT(P,TS,2)
         ELSE
            CALL FSATP(DV,DL,TOUT,TS,P)
            V1=1.D0/DL
            V2=1.D0/DV
            H1=HVT3N(V1,TS)
            H2=HVT3N(V2,TS)
         ENDIF
         X=(H-H1)/(H2-H1)
         VBPH=V1+X*(V2-V1)
      ELSE
         VBPH=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION SBPH(P,H)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF SPECIFIC VOLUME
C FOR GIVEN VALUES OF P AND H
C
C INPUT:   P  Pressure [MPa]
C          H     Specific Enthalpy [kilojoule per kilogram]
C
C OUTPUT:  S     Specific Entropy [kilojoule per kilogram and Kelvin]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSPH(P,H,IREG)
C  
      IF ((IREG .EQ. 1) .OR. (IREG .EQ. 2) .OR. 
     1    (IREG .EQ. 3) .OR. (IREG .EQ. 5)) THEN
         TB=TPH(P,H,IREG)
         SBPH=SPT(P,TB,IREG)
      ELSEIF (IREG .EQ. 9) THEN
         TS=TSATPN(P)
         IF (TS .LE. 623.15D0) THEN
            S1=SPT(P,TS,1)
            S2=SPT(P,TS,2)
            H1=HPT(P,TS,1)
            H2=HPT(P,TS,2)
         ELSE
            CALL FSATP(DV,DL,TOUT,TS,P)
            V1=1.D0/DL
            V2=1.D0/DV
            H1=HVT3N(V1,TS)
            H2=HVT3N(V2,TS)
            S1=SVT3N(V1,TS)
            S2=SVT3N(V2,TS)
         ENDIF
         X=(H-H1)/(H2-H1)
         SBPH=S1+X*(S2-S1)
      ELSE
         SBPH=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION CPBPH(P,H)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF       SPECIFIC ISOBARIC HEAT CAPACITY
C FOR GIVEN VALUES OF P AND H
C
C INPUT:   P  Pressure [MPa]
C          H     Specific Enthalpy [kilojoule per kilogram]
C
C OUTPUT:  CP    Specific Isobaric Heat Capacity 
C                [kilojoule per kilogram and Kelvin]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSPH(P,H,IREG)
C  
      IF ((IREG .EQ. 1) .OR. (IREG .EQ. 2) .OR. 
     1    (IREG .EQ. 3) .OR. (IREG .EQ. 5)) THEN
         T=TPH(P,H,IREG)
         CPBPH=CPPT(P,T,IREG)
      ELSE
         CPBPH=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION WBPH(P,H)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF SPEED OF SOUND
C FOR GIVEN VALUES OF P AND H
C
C INPUT:   P  Pressure [MPa]
C          H     Specific Enthalpy [kilojoule per kilogram]
C
C OUTPUT:  W     Speed of Sound [meter per second]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSPH(P,H,IREG)
C  
      IF ((IREG .EQ. 1) .OR. (IREG .EQ. 2) .OR. 
     1    (IREG .EQ. 3) .OR. (IREG .EQ. 5)) THEN
         T=TPH(P,H,IREG)
         WBPH=WPT(P,T,IREG)
      ELSE
         WBPH=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION XBPH(P,H)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF DRYNESS FRACTION
C FOR GIVEN VALUES OF P AND H
C
C INPUT:   P  Pressure [MPa]
C          H     Specific Enthalpy [kilojoule per kilogram]
C
C OUTPUT:  X     Dryness Fraction [-]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSPH(P,H,IREG)
C  
      IF (IREG .EQ. 9) THEN
         TS=TSATPN(P)
         IF (TS .LE. 623.15D0) THEN
            H1=HPT1N(P,TS)
            H2=HPT2N(P,TS)
         ELSE
            CALL FSATP(DV,DL,TOUT,TS,P)
            V1=1.D0/DL
            V2=1.D0/DV
            H1=HVT3N(V1,TS)
            H2=HVT3N(V2,TS)
         ENDIF
         XBPH=(H-H1)/(H2-H1)
      ELSE
         XBPH=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION TBPS(P,S)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF TEMPERATURE
C FOR GIVEN VALUES OF P AND S
C
C INPUT:   P  Pressure [MPa]
C          S     Specific Entropy [kilojoule per kilogram and Kelvin]
C
C OUTPUT:  T     Temperature [Kelvin]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSPS(P,S,IREG)
C  
      IF ((IREG .EQ. 1) .OR. (IREG .EQ. 2) .OR. 
     1    (IREG .EQ. 3) .OR. (IREG .EQ. 5)) THEN
         TBPS=TPS(P,S,IREG)
      ELSEIF (IREG .EQ. 9) THEN
         TBPS=TSATPN(P)
      ELSE
         TBPS=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION VBPS(P,S)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF SPECIFIC VOLUME
C FOR GIVEN VALUES OF P AND S
C
C INPUT:   P  Pressure [MPa]
C          S     Specific Entropy [kilojoule per kilogram and Kelvin]
C
C OUTPUT:  V     Specific Volume [cubic meter per kilogram]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSPS(P,S,IREG)
C
      IF ((IREG .EQ. 1) .OR. (IREG .EQ. 2) .OR. 
     1    (IREG .EQ. 3) .OR. (IREG .EQ. 5)) THEN
         TB=TPS(P,S,IREG)
         VBPS=VPT(P,TB,IREG)
      ELSEIF (IREG .EQ. 9) THEN
         TS=TSATPN(P)
         IF (TS .LE. 623.15D0) THEN
            V1=VPT(P,TS,1)
            V2=VPT(P,TS,2)
            S1=SPT(P,TS,1)
            S2=SPT(P,TS,2)
         ELSE
            CALL FSATP(DV,DL,TOUT,TS,P)
            V1=1.D0/DL
            V2=1.D0/DV
            S1=SVT3N(V1,TS)
            S2=SVT3N(V2,TS)
         ENDIF
         X=(S-S1)/(S2-S1)
         VBPS=V1+X*(V2-V1)
      ELSE
         VBPS=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION HBPS(P,S)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF SPECIFIC ENTHALPY
C FOR GIVEN VALUES OF P AND S
C
C INPUT:   P  Pressure [MPa]
C          S     Specific Entropy [kilojoule per kilogram and Kelvin]
C
C OUTPUT:  H     Specific Enthalpy [kilojoule per kilogram]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSPS(P,S,IREG)
C  
      IF ((IREG .EQ. 1) .OR. (IREG .EQ. 2) .OR. 
     1    (IREG .EQ. 3) .OR. (IREG .EQ. 5)) THEN
         TB=TPS(P,S,IREG)
         HBPS=HPT(P,TB,IREG)
      ELSEIF (IREG .EQ. 9) THEN
         TS=TSATPN(P)
         IF (TS .LE. 623.15D0) THEN
            S1=SPT(P,TS,1)
            S2=SPT(P,TS,2)
            H1=HPT(P,TS,1)
            H2=HPT(P,TS,2)
         ELSE
            CALL FSATP(DV,DL,TOUT,TS,P)
            V1=1.D0/DL
            V2=1.D0/DV
            H1=HVT3N(V1,TS)
            H2=HVT3N(V2,TS)
            S1=SVT3N(V1,TS)
            S2=SVT3N(V2,TS)
         ENDIF
         X=(S-S1)/(S2-S1)
         HBPS=H1+X*(H2-H1)
      ELSE
         HBPS=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION CPBPS(P,S)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF       SPECIFIC ISOBARIC HEAT CAPACITY
C FOR GIVEN VALUES OF P AND S
C
C INPUT:   P  Pressure [MPa]
C          S     Specific Entropy [kilojoule per kilogram and Kelvin]
C
C OUTPUT:  CP    Specific Isobaric Heat Capacity 
C                [kilojoule per kilogram and Kelvin]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSPS(P,S,IREG)
C  
      IF ((IREG .EQ. 1) .OR. (IREG .EQ. 2) .OR. 
     1    (IREG .EQ. 3) .OR. (IREG .EQ. 5)) THEN
         T=TPS(P,S,IREG)
         CPBPS=CPPT(P,T,IREG)
      ELSE
         CPBPS=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION WBPS(P,S)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF       SPEED OF SOUND
C FOR GIVEN VALUES OF P AND S
C
C INPUT:   P  Pressure [MPa]
C          S     Specific Entropy [kilojoule per kilogram and Kelvin]
C
C OUTPUT:  W     Speed of Sound [meter per second]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSPS(P,S,IREG)
C  
      IF ((IREG .EQ. 1) .OR. (IREG .EQ. 2) .OR. 
     1    (IREG .EQ. 3) .OR. (IREG .EQ. 5)) THEN
         T=TPS(P,S,IREG)
         WBPS=WPT(P,T,IREG)
      ELSE
         WBPS=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION XBPS(P,S)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF DRYNESS FRACTION
C FOR GIVEN VALUES OF P AND S
C
C INPUT:   P  Pressure [MPa]
C          S     Specific Entropy [kilojoule per kilogram and Kelvin]
C
C OUTPUT:  X     Dryness Fraction [-]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSPS(P,S,IREG)
C  
      IF (IREG .EQ. 9) THEN
         TS=TSATPN(P)
         IF (TS .LE. 623.15D0) THEN
            S1=SPT1N(P,TS)
            S2=SPT2N(P,TS)
         ELSE
            CALL FSATP(DV,DL,TOUT,TS,P)
            V1=1.D0/DL
            V2=1.D0/DV
            S1=SVT3N(V1,TS)
            S2=SVT3N(V2,TS)
         ENDIF
         XBPS=(S-S1)/(S2-S1)
      ELSE
         XBPS=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION TBPV(P,V)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF TEMPERATURE
C FOR GIVEN VALUES OF P AND V
C
C INPUT:   P  Pressure [MPa]
C          V     Specific Volume [cubic meter per kilogram]
C
C OUTPUT:  T     Temperature [Kelvin]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSPV(P,V,1,IREG)
C  
      IF ((IREG .EQ. 1) .OR. (IREG .EQ. 2) .OR. 
     1    (IREG .EQ. 3) .OR. (IREG .EQ. 5)) THEN
         CALL TPV(P,V,IREG,1,T1,T2)
         TBPV=T1
      ELSEIF (IREG .EQ. 9) THEN
         TBPV=TSATPN(P)
      ELSE
         TBPV=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION TBPVZ(P,V)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF TEMPERATURE
C FOR GIVEN VALUES OF P AND V (additional Function)
C
C INPUT:   P  Pressure [MPa]
C          V     Specific Volume [cubic meter per kilogram]
C
C OUTPUT:  T     Temperature [Kelvin]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
        TGR=TPVGR(P)
        IF ((TGR .GE. 273.15D0) .AND. (TGR .LE. 277.13D0)) THEN
C
                CALL REGSPV(P,V,2,IREG)
C  
                IF ((IREG .EQ. 1) .OR. (IREG .EQ. 2) .OR. 
     1          (IREG .EQ. 3) .OR. (IREG .EQ. 5)) THEN
                CALL TPV(P,V,IREG,2,T1,T2)
                        TBPVZ=T1
                ELSEIF (IREG .EQ. 9) THEN
                        TBPVZ=TSATPN(P)
                ELSE
                        TBPVZ=0.D0
                ENDIF
C
        ELSE
                TBPVZ=0.D0
        ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION HBPV(P,V)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF SPECIFIC ENTHALPY
C FOR GIVEN VALUES OF P AND V
C
C INPUT:   P  Pressure [MPa]
C          V     Specific Volume [cubic meter per kilogram]
C
C OUTPUT:  H     Specific Enthalpy [kilojoule per kilogram]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSPV(P,V,1,IREG)
C
      IF ((IREG .EQ. 1) .OR. (IREG .EQ. 2) .OR. 
     1    (IREG .EQ. 3) .OR. (IREG .EQ. 5)) THEN
         CALL TPV(P,V,IREG,1,T1,T2)
           TB=T1
         HBPV=HPT(P,TB,IREG)
      ELSEIF (IREG .EQ. 9) THEN
         TS=TSATPN(P)
         IF (TS .LE. 623.15D0) THEN
            V1=VPT(P,TS,1)
            V2=VPT(P,TS,2)
            H1=HPT(P,TS,1)
            H2=HPT(P,TS,2)
         ELSE
            CALL FSATP(DV,DL,TOUT,TS,P)
            V1=1.D0/DL
            V2=1.D0/DV
            H1=HVT3N(V1,TS)
            H2=HVT3N(V2,TS)
         ENDIF
         X=(V-V1)/(V2-V1)
         HBPV=H1+X*(H2-H1)
      ELSE
         HBPV=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION SBPV(P,V)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF SPECIFIC ENTROPY
C FOR GIVEN VALUES OF P AND V
C
C INPUT:   P  Pressure [MPa]
C          V     Specific Volume [cubic meter per kilogram]
C
C OUTPUT:  S     Specific Entropy [kilojoule per kilogram and Kelvin]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSPV(P,V,1,IREG)
C  
      IF ((IREG .EQ. 1) .OR. (IREG .EQ. 2) .OR. 
     1    (IREG .EQ. 3) .OR. (IREG .EQ. 5)) THEN
         CALL TPV(P,V,IREG,1,T1,T2)
           TB=T1
         SBPV=SPT(P,TB,IREG)
      ELSEIF (IREG .EQ. 9) THEN
         TS=TSATPN(P)
         IF (TS .LE. 623.15D0) THEN
            S1=SPT(P,TS,1)
            S2=SPT(P,TS,2)
            V1=VPT(P,TS,1)
            V2=VPT(P,TS,2)
         ELSE
            CALL FSATP(DV,DL,TOUT,TS,P)
            V1=1.D0/DL
            V2=1.D0/DV
            S1=SVT3N(V1,TS)
            S2=SVT3N(V2,TS)
         ENDIF
         X=(V-V1)/(V2-V1)
         SBPV=S1+X*(S2-S1)
      ELSE
         SBPV=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION CPBPV(P,V)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF       SPECIFIC ISOBARIC HEAT CAPACITY
C FOR GIVEN VALUES OF P AND V
C
C INPUT:   P  Pressure [MPa]
C          V     Specific Volume [cubic meter per kilogram]
C
C OUTPUT:  CP    Specific Isobaric Heat Capacity 
C                [kilojoule per kilogram and Kelvin]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSPV(P,V,1,IREG)
C  
      IF ((IREG .EQ. 1) .OR. (IREG .EQ. 2) .OR. 
     1    (IREG .EQ. 3) .OR. (IREG .EQ. 5)) THEN
         CALL TPV(P,V,IREG,1,T1,T2)
           T=T1
         CPBPV=CPPT(P,T,IREG)
      ELSE
         CPBPV=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION WBPV(P,V)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF       SPEED OF SOUND
C FOR GIVEN VALUES OF P AND V
C
C INPUT:   P  Pressure [MPa]
C          V     Specific Volume [cubic meter per kilogram]
C
C OUTPUT:  W     Speed of Sound [meter per second]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSPV(P,V,1,IREG)
C  
      IF ((IREG .EQ. 1) .OR. (IREG .EQ. 2) .OR. 
     1    (IREG .EQ. 3) .OR. (IREG .EQ. 5)) THEN
         CALL TPV(P,V,IREG,1,T1,T2)
           T=T1
         WBPV=WPT(P,T,IREG)
      ELSE
         WBPV=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION XBPV(P,V)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF DRYNESS FRACTION
C FOR GIVEN VALUES OF P AND V
C
C INPUT:   P  Pressure [MPa]
C          V     Specific Volume [cubic meter per kilogram]
C
C OUTPUT:  X     Dryness Fraction [-]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSPV(P,V,1,IREG)
C  
      IF (IREG .EQ. 9) THEN
         TS=TSATPN(P)
         IF (TS .LE. 623.15D0) THEN
            V1=VPT1N(P,TS)
            V2=VPT2N(P,TS)
         ELSE
            CALL FSATP(DV,DL,TOUT,TS,P)
            V1=1.D0/DL
            V2=1.D0/DV
         ENDIF
         XBPV=(V-V1)/(V2-V1)
      ELSE
         XBPV=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION PBTH(T,H)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF PRESSURE
C FOR GIVEN VALUES OF T AND H
C
C INPUT:   T  Temperature [Kelvin]
C          H     Specific Enthalpy [kilojoule per kilogram]
C
C OUTPUT:  P     Pressure [MPa]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSTH(T,H,1,IREG)
C  
      IF ((IREG .EQ. 1) .OR. (IREG .EQ. 2) .OR. 
     1    (IREG .EQ. 3) .OR. (IREG .EQ. 5)) THEN
         CALL PTH(T,H,IREG,1,P1,P2)
         PBTH=P1
      ELSEIF (IREG .EQ. 9) THEN
         PBTH=PSATTN(T)
      ELSE
         PBTH=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION PBTHZ(T,H)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF PRESSURE
C FOR GIVEN VALUES OF T AND H (additional Function)
C
C INPUT:   T  Temperature [Kelvin]
C          H     Specific Enthalpy [kilojoule per kilogram]
C
C OUTPUT:  P     Pressure [MPa]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      IF (T .LE. 611.9697430625D0) THEN
                CALL REGSTH(T,H,2,IREG)
C  
                IF ((IREG .EQ. 1) .OR. (IREG .EQ. 2) .OR. 
     1          (IREG .EQ. 3) .OR. (IREG .EQ. 5)) THEN
                        CALL PTH(T,H,IREG,2,P1,P2)
                        PBTHZ=P1
                ELSEIF (IREG .EQ. 9) THEN
                        PBTHZ=PSATTN(T)
                ELSE
                        PBTHZ=0.D0
                ENDIF
C
        ELSE
                PBTHZ=0.D0
        ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION PBTS(T,S)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF PRESSURE
C FOR GIVEN VALUES OF T AND S
C
C INPUT:   T  Temperature [Kelvin]
C          S     Specific Entropy [kilojoule per kilogram and Kelvin]
C
C OUTPUT:  P     Pressure [MPa]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSTS(T,S,1,IREG)
C  
      IF ((IREG .EQ. 1) .OR. (IREG .EQ. 2) .OR. 
     1    (IREG .EQ. 3) .OR. (IREG .EQ. 5)) THEN
         CALL PTS(T,S,IREG,1,P1,P2)
         PBTS=P1
      ELSEIF (IREG .EQ. 9) THEN
         PBTS=PSATTN(T)
      ELSE
         PBTS=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION PBTSZ(T,S)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF PRESSURE
C FOR GIVEN VALUES OF T AND S (additional Function)
C
C INPUT:   T  Temperature [Kelvin]
C          S     Specific Entropy [kilojoule per kilogram and Kelvin]
C
C OUTPUT:  P     Pressure [MPa]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      IF ((T .GE. 273.15D0) .AND. (T .LE. 277.13D0)) THEN
                CALL REGSTS(T,S,2,IREG)
C
                IF ((IREG .EQ. 1) .OR. (IREG .EQ. 2) .OR. 
     1          (IREG .EQ. 3) .OR. (IREG .EQ. 5)) THEN
                        CALL PTS(T,S,IREG,2,P1,P2)
                        PBTSZ=P1
                ELSEIF (IREG .EQ. 9) THEN
                        PBTSZ=PSATTN(T)
                ELSE
                        PBTSZ=0.D0
                ENDIF
C
        ELSE
                PBTSZ=0.D0
        ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION PBTV(T,V)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF PRESSURE
C FOR GIVEN VALUES OF T AND V
C
C INPUT:   T  Temperature [Kelvin]
C          V     Specific Volume [cubic meter per kilogram]
C
C OUTPUT:  P     Pressure [MPa]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSTV(T,V,IREG)
C  
      IF ((IREG .EQ. 1) .OR. (IREG .EQ. 2) .OR. 
     1    (IREG .EQ. 3) .OR. (IREG .EQ. 5)) THEN
         PBTV=PTV(T,V,IREG)
      ELSEIF (IREG .EQ. 9) THEN
         PBTV=PSATTN(T)
      ELSE
         PBTV=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION PBVH(V,H)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF PRESSURE
C FOR GIVEN VALUES OF V AND H
C
C INPUT:   V     Specific Volume [cubic meter per kilogram]
C          H     Specific Enthalpy [kilojoule per kilogram]
C
C OUTPUT:  P     Pressure [MPa]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSHV(H,V,IREG)
C  
      IF (IREG .NE. 0)  THEN
         PBVH=PHV(H,V,IREG)
      ELSE
         PBVH=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION PBVS(V,S)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF PRESSURE
C FOR GIVEN VALUES OF V AND S
C
C INPUT:   V     Specific Volume [cubic meter per kilogram]
C          S     Specific Entropy [kilojoule per kilogram and Kelvin]
C
C OUTPUT:  P     Pressure [MPa]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSSV(S,V,IREG)
C  
      IF (IREG .NE. 0)  THEN
         PBVS=PSV(S,V,IREG)
      ELSE
         PBVS=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION PBHS(H,S)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF PRESSURE
C FOR GIVEN VALUES OF H AND S
C
C INPUT:   H     Specific Enthalpy [kilojoule per kilogram]
C          S     Specific Entropy [kilojoule per kilogram and Kelvin]
C
C OUTPUT:  P     Pressure [MPa]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSHS(H,S,IREG)
C  
      IF (IREG .NE. 0)  THEN
         PBHS=PHS(H,S,IREG)
      ELSE
         PBHS=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION TBHS(H,S)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF TEMPERATURE
C FOR GIVEN VALUES OF H AND S
C
C INPUT:   H     Specific Enthalpy [kilojoule per kilogram]
C          S     Specific Entropy [kilojoule per kilogram and Kelvin]
C
C OUTPUT:  P     Temperature [Kelvin]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      CALL REGSHS(H,S,IREG)
C  
      IF (IREG .NE. 0)  THEN
         P=PHS(H,S,IREG)
           TBHS=TPH(P,H,IREG)
      ELSE
         TBHS=0.D0
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION VBPX(P,X)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF SPECIFIC VOLUME
C FOR GIVEN VALUES OF P AND X
C
C INPUT:   P  Pressure [MPa]
C          X     Dryness Fraction [-]
C
C OUTPUT:  V     Specific Volume [cubic meter per kilogram]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      IF ((X .GT. 1.D0) .OR. (X .LT. 0.D0) .OR. (P .GT. 22.064D0)) THEN
         VBPX=0.D0
      ELSE
         TS=TSATPN(P)
         IF (TS .LE. 623.15D0) THEN
            V1=VPT(P,TS,1)
            V2=VPT(P,TS,2)
         ELSE
            CALL FSATP(DV,DL,TOUT,TS,P)
            V1=1.D0/DL
            V2=1.D0/DV
         ENDIF
         VBPX=V1+X*(V2-V1)
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION VBTX(T,X)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF SPECIFIC VOLUME
C FOR GIVEN VALUES OF T AND X
C
C INPUT:   T  Temperature [Kelvin]
C          X     Dryness Fraction [-]
C
C OUTPUT:  V     Specific Volume [cubic meter per kilogram]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      IF ((X .GT. 1.D0) .OR. (X .LT. 0.D0) .OR. 
     &(T .LT. 273.15D0) .OR. (T .GT. 647.096D0)) THEN
         VBTX=0.D0
      ELSE
         PS=PSATTN(T)
         IF (T .LE. 623.15D0) THEN
            V1=VPT(PS,T,1)
            V2=VPT(PS,T,2)
         ELSE
            CALL FSATP(DV,DL,TOUT,T,PS)
            V1=1.D0/DL
            V2=1.D0/DV
         ENDIF
         VBTX=V1+X*(V2-V1)
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION HBPX(P,X)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF SPECIFIC ENTHALPY
C FOR GIVEN VALUES OF P AND X
C
C INPUT:   P  Pressure [MPa]
C          X     Dryness Fraction [-]
C
C OUTPUT:  H     Specific Enthalpy [kilojoule per kilogram]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      IF ((X .GT. 1.D0) .OR. (X .LT. 0.D0) .OR. (P .GT. 22.064D0)) THEN
         HBPX=0.D0
      ELSE
         TS=TSATPN(P)
         IF (TS .LE. 623.15D0) THEN
            H1=HPT(P,TS,1)
            H2=HPT(P,TS,2)
         ELSE
            CALL FSATP(DV,DL,TOUT,TS,P)
            V1=1.D0/DL
            V2=1.D0/DV
            H1=HVT3N(V1,TS)
            H2=HVT3N(V2,TS)
         ENDIF
         HBPX=H1+X*(H2-H1)
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION HBTX(T,X)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF SPECIFIC ENTHALPY
C FOR GIVEN VALUES OF T AND X
C
C INPUT:   T  Temperature [Kelvin]
C          X     Dryness Fraction [-]
C
C OUTPUT:  H     Specific Enthalpy [kilojoule per kilogram]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      IF ((X .GT. 1.D0) .OR. (X .LT. 0.D0) .OR. 
     &(T .LT. 273.15D0) .OR. (T .GT. 647.096D0)) THEN
         HBTX=0.D0
      ELSE
         PS=PSATTN(T)
         IF (T .LE. 623.15D0) THEN
            H1=HPT(PS,T,1)
            H2=HPT(PS,T,2)
         ELSE
            CALL FSATP(DV,DL,TOUT,T,PS)
            V1=1.D0/DL
            V2=1.D0/DV
            H1=HVT3N(V1,T)
            H2=HVT3N(V2,T)
         ENDIF
         HBTX=H1+X*(H2-H1)
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION SBPX(P,X)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF SPECIFIC ENTRPOPY
C FOR GIVEN VALUES OF P AND X
C
C INPUT:   P  Pressure [MPa]
C          X     Dryness Fraction [-]
C
C OUTPUT:  S     Specific Entropy [kilojoule per kilogram and Kelvin]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      IF ((X .GT. 1.D0) .OR. (X .LT. 0.D0) .OR. (P .GT. 22.064D0)) THEN
         SBPX=0.D0
      ELSE
         TS=TSATPN(P)
         IF (TS .LE. 623.15D0) THEN
            S1=SPT(P,TS,1)
            S2=SPT(P,TS,2)
         ELSE
            CALL FSATP(DV,DL,TOUT,TS,P)
            V1=1.D0/DL
            V2=1.D0/DV
            S1=SVT3N(V1,TS)
            S2=SVT3N(V2,TS)
         ENDIF
         SBPX=S1+X*(S2-S1)
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION SBTX(T,X)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF SPECIFIC ENTRPOPY
C FOR GIVEN VALUES OF T AND X
C
C INPUT:   T  Temperature [Kelvin]
C          X     Dryness Fraction [-]
C
C OUTPUT:  S     Specific Entropy [kilojoule per kilogram and Kelvin]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      IF ((X .GT. 1.D0) .OR. (X .LT. 0.D0) .OR. 
     &(T .LT. 273.15D0) .OR. (T .GT. 647.096D0)) THEN
         SBTX=0.D0
      ELSE
         PS=PSATTN(T)
         IF (T .LE. 623.15D0) THEN
            S1=SPT(PS,T,1)
            S2=SPT(PS,T,2)
         ELSE
            CALL FSATP(DV,DL,TOUT,T,PS)
            V1=1.D0/DL
            V2=1.D0/DV
            S1=SVT3N(V1,T)
            S2=SVT3N(V2,T)
         ENDIF
         SBTX=S1+X*(S2-S1)
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION PBT(T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF       SATURATION PRESSURE
C FOR A GIVEN VALUE OF T
C
C INPUT:   T  Temperature [Kelvin]
C
C OUTPUT:  P     Saturation Pressure [MPa]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      IF ((T .LE. 647.096D0) .AND. (T .GE. 273.15D0)) THEN
         PBT=PSATTN(T)
      ELSE
         PBT=0.D0
      ENDIF   
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION TBP(P)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF       SATURATION TEMPERATURE
C FOR A GIVEN VALUE OF P
C
C INPUT:   P  Pressure [MPa]
C
C OUTPUT:  T     Temperature [Kelvin]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      IF ((P .LE. 22.064D0) .AND. (P .GE. 5.D-4)) THEN
         TBP=TSATPN(P)
      ELSE
         TBP=0.D0
      ENDIF   
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION H1BP(P)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF       SPECIFIC ENTHALPY ON THE SATURATED 
C LIQUID LINE
C FOR A GIVEN VALUE OF P
C
C INPUT:   P  Pressure [MPa]
C
C OUTPUT:  H     Specific Enthalpy on the Saturated Liquid Line
C                [kilojoule per kilogram]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      IF ((P .LE. 22.064D0) .AND. (P .GE. 5.D-4)) THEN
         TS=TSATPN(P)
         IF (TS .LE. 623.15D0) THEN
            H1BP=HPT1N(P,TS)
         ELSE
            CALL FSATP(DV,DL,TOUT,TS,P)
            V1=1.D0/DL
            H1BP=HVT3N(V1,TS)
         ENDIF
      ELSE
         H1BP=0.D0
      ENDIF   
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION H2BP(P)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF       SPECIFIC ENTHALPY ON THE SATURATED 
C VAPOUR LINE
C FOR A GIVEN VALUE OF P
C
C INPUT:   P  Pressure [MPa]
C
C OUTPUT:  H     Specific Enthalpy on the Saturated Vapour Line 
C                [kilojoule per kilogram]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      IF ((P .LE. 22.064D0) .AND. (P .GE. 5.D-4)) THEN
         TS=TSATPN(P)
         IF (TS .LE. 623.15D0) THEN
            H2BP=HPT2N(P,TS)
         ELSE
            CALL FSATP(DV,DL,TOUT,TS,P)
            V2=1.D0/DV
            H2BP=HVT3N(V2,TS)
         ENDIF
      ELSE
         H2BP=0.D0
      ENDIF   
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION S1BP(P)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF       SPECIFIC ENTROPY ON THE SATURATED 
C LIQUID LINE
C FOR A GIVEN VALUE OF P
C
C INPUT:   P  Pressure [MPa]
C
C OUTPUT:  S     Specific Entropy on the Saturated Liquid Line
C                [kilojoule per kilogram and Kelvin]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      IF ((P .LE. 22.064D0) .AND. (P .GE. 5.D-4)) THEN
         TS=TSATPN(P)
         IF (TS .LE. 623.15D0) THEN
            S1BP=SPT1N(P,TS)
         ELSE
            CALL FSATP(DV,DL,TOUT,TS,P)
            V1=1.D0/DL
            S1BP=SVT3N(V1,TS)
         ENDIF
      ELSE
         S1BP=0.D0
      ENDIF   
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION S2BP(P)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF       SPECIFIC ENTROPY ON THE SATURATED 
C VAPOUR LINE
C FOR A GIVEN VALUE OF P
C
C INPUT:   P  Pressure [MPa]
C
C OUTPUT:  S     Specific Entropy on the Saturated Vapour Line 
C                [kilojoule per kilogram and Kelvin]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      IF ((P .LE. 22.064D0) .AND. (P .GE. 5.D-4)) THEN
         TS=TSATPN(P)
         IF (TS .LE. 623.15D0) THEN
            S2BP=SPT2N(P,TS)
         ELSE
            CALL FSATP(DV,DL,TOUT,TS,P)
            V2=1.D0/DV
            S2BP=SVT3N(V2,TS)
         ENDIF
      ELSE
         S2BP=0.D0
      ENDIF   
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION V1BP(P)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF       SPECIFIC VOLUME ON THE SATURATED 
C LIQUID LINE
C FOR A GIVEN VALUE OF P
C
C INPUT:   P  Pressure [MPa]
C
C OUTPUT:  V     Specific Volume on the Saturated Liquid Line
C                [cubic meter per kilogram]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      IF ((P .LE. 22.064D0) .AND. (P .GE. 5.D-4)) THEN
         TS=TSATPN(P)
         IF (TS .LE. 623.15D0) THEN
            V1BP=VPT1N(P,TS)
         ELSE
            CALL FSATP(DV,DL,TOUT,TS,P)
            V1BP=1.D0/DL
         ENDIF
      ELSE
         V1BP=0.D0
      ENDIF   
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION V2BP(P)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF       SPECIFIC VOLUME ON THE SATURATED 
C VAPOUR LINE
C FOR A GIVEN VALUE OF P
C
C INPUT:   P  Pressure [MPa]
C
C OUTPUT:  V     Specific Volume on the Saturated Vapour Line 
C                [cubic meter per kilogram]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      IF ((P .LE. 22.064D0) .AND. (P .GE. 5.D-4)) THEN
         TS=TSATPN(P)
         IF (TS .LE. 623.15D0) THEN
            V2BP=VPT2N(P,TS)
         ELSE
            CALL FSATP(DV,DL,TOUT,TS,P)
            V2BP=1.D0/DV
         ENDIF
      ELSE
         V2BP=0.D0
      ENDIF   
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION CP1BP(P)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF       SPECIFIC ISOBARIC HEAT CAPACITY ON THE 
C SATURATED LIQUID LINE
C FOR A GIVEN VALUE OF P
C
C INPUT:   P  Pressure [MPa]
C
C OUTPUT:  CP    Specific Isobaric Heat Capacity on the Saturated 
C                Liquid Line 
C                [kilojoule per kilogram and Kelvin]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      IF ((P .LE. 22.064D0) .AND. (P .GE. 5.D-4)) THEN
         TS=TSATPN(P)
         IF (TS .LE. 623.15D0) THEN
            CP1BP=CPPT1N(P,TS)
         ELSE
            CALL FSATP(DV,DL,TOUT,TS,P)
            V1=1.D0/DL
            CP1BP=CPVT3N(V1,TS)
         ENDIF
      ELSE
         CP1BP=0.D0
      ENDIF   
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION CP2BP(P)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF       SPECIFIC ISOBARIC HEAT CAPACITY ON 
C THE SATURATED VAPOUR LINE
C FOR A GIVEN VALUE OF P
C
C INPUT:   P  Pressure [MPa]
C
C OUTPUT:  CP    Specific Isobaric Heat Capacity on the Saturated 
C                Vapour Line 
C                [kilojoule per kilogram and Kelvin]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      IF ((P .LE. 22.064D0) .AND. (P .GE. 5.D-4)) THEN
         TS=TSATPN(P)
         IF (TS .LE. 623.15D0) THEN
            CP2BP=CPPT2N(P,TS)
         ELSE
            CALL FSATP(DV,DL,TOUT,TS,P)
            V2=1.D0/DV
            CP2BP=CPVT3N(V2,TS)
         ENDIF
      ELSE
         CP2BP=0.D0
      ENDIF   
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION W1BP(P)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF       SPEED OF SOUND ON THE SATURATED 
C LIQUID LINE
C FOR A GIVEN VALUE OF P
C
C INPUT:   P  Pressure [MPa]
C
C OUTPUT:  W     Speed of Sound on the Saturated Liquid Line
C                [meter per second]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      IF ((P .LE. 22.064D0) .AND. (P .GE. 5.D-4)) THEN
         TS=TSATPN(P)
         IF (TS .LE. 623.15D0) THEN
            W1BP=WPT1N(P,TS)
         ELSE
            CALL FSATP(DV,DL,TOUT,TS,P)
            V1=1.D0/DL
            W1BP=WVT3N(V1,TS)
         ENDIF
      ELSE
         W1BP=0.D0
      ENDIF   
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION W2BP(P)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF       SPEED OF SOUND ON THE SATURATED 
C VAPOUR LINE
C FOR A GIVEN VALUE OF P
C
C INPUT:   P  Pressure [MPa]
C
C OUTPUT:  W     Speed of Sound on the Saturated Vapour Line 
C                [meter per second]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      IF ((P .LE. 22.064D0) .AND. (P .GE. 5.D-4)) THEN
         TS=TSATPN(P)
         IF (TS .LE. 623.15D0) THEN
            W2BP=WPT2N(P,TS)
         ELSE
            CALL FSATP(DV,DL,TOUT,TS,P)
            V2=1.D0/DV
            W2BP=WVT3N(V2,TS)
         ENDIF
      ELSE
         W2BP=0.D0
      ENDIF   
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION H1BT(T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF       SPECIFIC ENTHALPY ON THE SATURATED 
C LIQUID LINE
C FOR A GIVEN VALUE OF T
C
C INPUT:   T  Temperature [Kelvin]
C
C OUTPUT:  H     Specific Enthalpy on the Saturated Liquid Line
C                [kilojoule per kilogram]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      IF ((T .LE. 647.096D0) .AND. (T .GE. 273.15D0)) THEN
         PS=PSATTN(T)
         IF (T .LE. 623.15D0) THEN
            H1BT=HPT1N(PS,T)
         ELSE
            CALL FSATP(DV,DL,TOUT,T,PS)
            V1=1.D0/DL
            H1BT=HVT3N(V1,T)
         ENDIF
      ELSE
         H1BT=0.D0
      ENDIF   
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION H2BT(T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF       SPECIFIC ENTHALPY ON THE SATURATED 
C VAPOUR LINE 
C FOR A GIVEN VALUE OF T
C
C INPUT:   T  Temperature [Kelvin]
C
C OUTPUT:  H     Specific Enthalpy on the Saturated Vapour Line 
C                [kilojoule per kilogram]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      IF ((T .LE. 647.096D0) .AND. (T .GE. 273.15D0)) THEN
         PS=PSATTN(T)
         IF (T .LE. 623.15D0) THEN
            H2BT=HPT2N(PS,T)
         ELSE
            CALL FSATP(DV,DL,TOUT,T,PS)
            V2=1.D0/DV
            H2BT=HVT3N(V2,T)
         ENDIF
      ELSE
         H2BT=0.D0
      ENDIF   
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION S1BT(T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF       SPECIFIC ENTROPY ON THE SATURATED 
C LIQUID LINE
C FOR A GIVEN VALUE OF T
C
C INPUT:   T  Temperature [Kelvin]
C
C OUTPUT:  S     Specific Entropy on the Saturated Liquid Line
C                [kilojoule per kilogram and Kelvin]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      IF ((T .LE. 647.096D0) .AND. (T .GE. 273.15D0)) THEN
         PS=PSATTN(T)
         IF (T .LE. 623.15D0) THEN
            S1BT=SPT1N(PS,T)
         ELSE
            CALL FSATP(DV,DL,TOUT,T,PS)
            V1=1.D0/DL
            S1BT=SVT3N(V1,T)
         ENDIF
      ELSE
         S1BT=0.D0
      ENDIF   
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION S2BT(T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF       SPECIFIC ENTROPY ON THE SATURATED 
C VAPOUR LINE 
C FOR A GIVEN VALUE OF T
C
C INPUT:   T  Temperature [Kelvin]
C
C OUTPUT:  S     Specific Entropy on the Saturated Vapour Line 
C                [kilojoule per kilogram and Kelvin]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      IF ((T .LE. 647.096D0) .AND. (T .GE. 273.15D0)) THEN
         PS=PSATTN(T)
         IF (T .LE. 623.15D0) THEN
            S2BT=SPT2N(PS,T)
         ELSE
            CALL FSATP(DV,DL,TOUT,T,PS)
            V2=1.D0/DV
            S2BT=SVT3N(V2,T)
         ENDIF
      ELSE
         S2BT=0.D0
      ENDIF   
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION V1BT(T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF       SPECIFIC VOLUME ON THE SATURATED 
C LIQUID LINE
C FOR A GIVEN VALUE OF T
C
C INPUT:   T  Temperature [Kelvin]
C
C OUTPUT:  V     Specific Volume on the Saturated Liquid Line
C                [cubic meter per kilogram]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      IF ((T .LE. 647.096D0) .AND. (T .GE. 273.15D0)) THEN
         PS=PSATTN(T)
         IF (T .LE. 623.15D0) THEN
            V1BT=VPT1N(PS,T)
         ELSE
            CALL FSATP(DV,DL,TOUT,T,PS)
            V1BT=1.D0/DL
         ENDIF
      ELSE
         V1BT=0.D0
      ENDIF   
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION V2BT(T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF       SPECIFIC VOLUME ON THE SATURATED 
C VAPOUR LINE 
C FOR A GIVEN VALUE OF T
C
C INPUT:   T  Temperature [Kelvin]
C
C OUTPUT:  V     Specific Volume on the Saturated Vapour Line 
C                [cubic meter per kilogram]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      IF ((T .LE. 647.096D0) .AND. (T .GE. 273.15D0)) THEN
         PS=PSATTN(T)
         IF (T .LE. 623.15D0) THEN
            V2BT=VPT2N(PS,T)
         ELSE
            CALL FSATP(DV,DL,TOUT,T,PS)
            V2BT=1.D0/DV
         ENDIF
      ELSE
         V2BT=0.D0
      ENDIF   
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION CP1BT(T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF       SPECIFIC ISOBARIC HEAT CAPACITY ON THE
C SATURATED LIQUID LINE 
C FOR A GIVEN VALUE OF T
C
C INPUT:   T  Temperature [Kelvin]
C
C OUTPUT:  CP    Specific Isobaric Heat Capacity on the Saturated 
C                Liquid Line [kilojoule per kilogram and Kelvin]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      IF ((T .LE. 647.096D0) .AND. (T .GE. 273.15D0)) THEN
         PS=PSATTN(T)
         IF (T .LE. 623.15D0) THEN
            CP1BT=CPPT1N(PS,T)
         ELSE
            CALL FSATP(DV,DL,TOUT,T,PS)
            V1=1.D0/DL
            CP1BT=CPVT3N(V1,T)
         ENDIF
      ELSE
         CP1BT=0.D0
      ENDIF   
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION CP2BT(T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF       SPECIFIC ISOBARIC HEAT CAPACITY ON THE
C SATURATED VAPOUR LINE 
C FOR A GIVEN VALUE OF T
C
C INPUT:   T  Temperature [Kelvin]
C
C OUTPUT:  CP    Specific Isobaric Heat Capacity on the Saturated 
C                Vapour Line [kilojoule per kilogram and Kelvin]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      IF ((T .LE. 647.096D0) .AND. (T .GE. 273.15D0)) THEN
         PS=PSATTN(T)
         IF (T .LE. 623.15D0) THEN
            CP2BT=CPPT2N(PS,T)
         ELSE
            CALL FSATP(DV,DL,TOUT,T,PS)
            V2=1.D0/DV
            CP2BT=CPVT3N(V2,T)
         ENDIF
      ELSE
         CP2BT=0.D0
      ENDIF   
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION W1BT(T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF       SPEED OF SOUND ON THE SATURATED 
C LIQUID LINE
C FOR A GIVEN VALUE OF T
C
C INPUT:   T  Temperature [Kelvin]
C
C OUTPUT:  W     Speed of Sound on the Saturated Liquid Line
C                [meter per second]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      IF ((T .LE. 647.096D0) .AND. (T .GE. 273.15D0)) THEN
         PS=PSATTN(T)
         IF (T .LE. 623.15D0) THEN
            W1BT=WPT1N(PS,T)
         ELSE
            CALL FSATP(DV,DL,TOUT,T,PS)
            V1=1.D0/DL
            W1BT=WVT3N(V1,T)
         ENDIF
      ELSE
         W1BT=0.D0
      ENDIF   
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION W2BT(T)
C***********************************************************************
C
C FUNCTION FOR THE CALCULATION OF       SPEED OF SOUND ON THE SATURATED 
C VAPOUR LINE 
C FOR A GIVEN VALUE OF T
C
C INPUT:   T  Temperature [Kelvin]
C
C OUTPUT:  W     Speed of Sound on the Saturated Vapour Line 
C                [meter per second]
C
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      IF ((T .LE. 647.096D0) .AND. (T .GE. 273.15D0)) THEN
         PS=PSATTN(T)
         IF (T .LE. 623.15D0) THEN
            W2BT=WPT2N(PS,T)
         ELSE
            CALL FSATP(DV,DL,TOUT,T,PS)
            V2=1.D0/DV
            W2BT=WVT3N(V2,T)
         ENDIF
      ELSE
         W2BT=0.D0
      ENDIF   
C
      END
C
C***********************************************************************
      SUBROUTINE REGSOPT(P,T,IREG)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      double precision TOLD, POLD
      SAVE TOLD, POLD
        INTEGER IREGOLD
        SAVE IREGOLD
      data TOLD,POLD,IREGOLD / -1.0d0, -1.0d0, -1 /
C
      IF ((DABS(T-TOLD).LT. 1.D-6) .AND.
     *   (DABS(P-POLD).LT. 1.D-6)) THEN
         IREG=IREGOLD
         GOTO 999
      END IF
C
      IF ((T .LT. 273.15D0) .OR. (T .GT. 2273.15D0)
     1   .OR. (P .GT. 100.D0) .OR. (P .LT. 5.D-4)) THEN
         IREG=0
         GOTO 999
      ENDIF
C
        IREG=0
C
      IF (T .LE. 623.15D0) THEN
         PG=PSATTN(T)
         IF (P .GE. PG) THEN
            IREG=1
         ELSE
            IREG=2
         ENDIF
      ELSE IF (T .LE. 863.15D0) THEN
         PG=FB23(T)
         IF (P .GT. PG) THEN
            IREG=3
         ELSE
            IREG=2
         ENDIF
      ELSE IF (T .LE. 1073.15D0) THEN
         IREG=2
      ELSE
         IF (P .LE. 10.D0) THEN
            IREG=5
         ELSE
            IREG=0
         ENDIF
      ENDIF
C             
999   CONTINUE
C      
      TOLD=T
      POLD=P
      IREGOLD=IREG
C
      END      
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION VPT(P,T,IREG)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      IF (IREG .EQ. 0) THEN
         CALL REGSOPT(P,T,IREG)
      ENDIF
C  
      IF (IREG .EQ. 1) THEN
         VPT=VPT1N(P,T)
      ELSE IF (IREG .EQ. 2) THEN
         VPT=VPT2N(P,T)
      ELSE IF (IREG .EQ. 3) THEN
         VPT=VPT3N(P,T)
      ELSE IF (IREG .EQ. 5) THEN
         VPT=VPT5N(P,T)
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION HPT(P,T,IREG)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      IF (IREG .EQ. 0) THEN
         CALL REGSOPT(P,T,IREG)
      ENDIF
C  
      IF (IREG .EQ. 1) THEN
         HPT=HPT1N(P,T)
      ELSE IF (IREG .EQ. 2) THEN
         HPT=HPT2N(P,T)
      ELSE IF (IREG .EQ. 3) THEN
         V=VPT(P,T,3)
         HPT=HVT3N(V,T)
      ELSE IF (IREG .EQ. 5) THEN
         HPT=HPT5N(P,T)
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION SPT(P,T,IREG)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      IF (IREG .EQ. 0) THEN
         CALL REGSOPT(P,T,IREG)
      ENDIF
C  
      IF (IREG .EQ. 1) THEN
         SPT=SPT1N(P,T)
      ELSE IF (IREG .EQ. 2) THEN
         SPT=SPT2N(P,T)
      ELSE IF (IREG .EQ. 3) THEN
         V=VPT(P,T,3)
         SPT=SVT3N(V,T)
      ELSE IF (IREG .EQ. 5) THEN
         SPT=SPT5N(P,T)
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION WPT(P,T,IREG)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      IF (IREG .EQ. 0) THEN
         CALL REGSOPT(P,T,IREG)
      ENDIF
C  
      IF (IREG .EQ. 1) THEN
         WPT=WPT1N(P,T)
      ELSE IF (IREG .EQ. 2) THEN
         WPT=WPT2N(P,T)
      ELSE IF (IREG .EQ. 3) THEN
         V=VPT(P,T,3)
         WPT=WVT3N(V,T)
      ELSE IF (IREG .EQ. 5) THEN
         WPT=WPT5N(P,T)
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION CPPT(P,T,IREG)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      IF (IREG .EQ. 0) THEN
         CALL REGSOPT(P,T,IREG)
      ENDIF
C  
      IF (IREG .EQ. 1) THEN
         CPPT=CPPT1N(P,T)
      ELSE IF (IREG .EQ. 2) THEN
         CPPT=CPPT2N(P,T)
      ELSE IF (IREG .EQ. 3) THEN
         V=VPT(P,T,3)
         CPPT=CPVT3N(V,T)
      ELSE IF (IREG .EQ. 5) THEN
         CPPT=CPPT5N(P,T)
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION CVPT(P,T,IREG)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      IF (IREG .EQ. 0) THEN
         CALL REGSOPT(P,T,IREG)
      ENDIF
C  
      IF (IREG .EQ. 1) THEN
         CVPT=CVPT1N(P,T)
      ELSE IF (IREG .EQ. 2) THEN
         CVPT=CVPT2N(P,T)
      ELSE IF (IREG .EQ. 3) THEN
         V=VPT(P,T,3)
         CVPT=CVVT3N(V,T)
      ELSE IF (IREG .EQ. 5) THEN
         CVPT=CVPT5N(P,T)
      ENDIF
C
      END
C***********************************************************************
      DOUBLE PRECISION FUNCTION UPT(P,T,IREG)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      IF (IREG .EQ. 0) THEN
         CALL REGSOPT(P,T,IREG)
      ENDIF
C  
      IF (IREG .EQ. 1) THEN
         UPT=UPT1N(P,T)
      ELSE IF (IREG .EQ. 2) THEN
         UPT=UPT2N(P,T)
      ELSE IF (IREG .EQ. 3) THEN
         V=VPT3N(P,T)
         UPT=UVT3N(V,T)
      ELSE IF (IREG .EQ. 5) THEN
         UPT=UPT5N(P,T)
      ELSE
         UPT=0.D0
      ENDIF
C
      END
C***********************************************************************
      DOUBLE PRECISION FUNCTION GPT(P,T,IREG)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      IF (IREG .EQ. 0) THEN
         CALL REGSOPT(P,T,IREG)
      ENDIF
C  
      IF (IREG .EQ. 1) THEN
         GPT=GPT1N(P,T)
      ELSE IF (IREG .EQ. 2) THEN
         GPT=GPT2N(P,T)
      ELSE IF (IREG .EQ. 3) THEN
         V=VPT3N(P,T)
         GPT=GVT3N(V,T)
      ELSE IF (IREG .EQ. 5) THEN
         GPT=GPT5N(P,T)
      ELSE
         GPT=0.D0
      ENDIF
C
      END
C************************************************************************
      DOUBLEPRECISION FUNCTION VPT3N(P,T)
C************************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      double precision POLD, TOLD, VPT3NOLD
      SAVE POLD, TOLD, VPT3NOLD
      data POLD, TOLD, VPT3NOLD / -1.D+0, -1.D+0, -1.D+0 /
C
      IF ((DABS(P-POLD).LT. 1.D-6) .AND.
     *   (DABS(T-TOLD).LT. 1.D-6)) THEN
         VPT3N=VPT3NOLD
         GOTO 999
      END IF
C
      EPS=1.D-6
C
      CALL REG3S(P,T,IREG3)
C
      DEST=1.D0/VEST3(P,T,IREG3)
C
      VPT3N=1.D0/DITER3(P,T,DEST,EPS)
C
999     CONTINUE
C
      POLD=P
      TOLD=T
      VPT3NOLD=VPT3N
C
      END
C
C***********************************************************************
      SUBROUTINE REG3S(P,T,IREG3)
C***********************************************************************
        IMPLICIT double precision(A-H,O-Z)
C
      COMMON/CSUB2/TC,PC,DC
C
      IF (P .LE. 40.D0) THEN
C      
         IF (P .LE. PC) THEN
            TS = TSATPN(P)
            IF (T .LT. TS) THEN
               IREG3=1
            ELSE 
               IREG3=4
            ENDIF
         ELSE
            PG=PC+(T-TC)/3.727888D0
            IF (P .LE. 24.D0) THEN
               IF (T. LE. TC) THEN
                  IREG3=1
               ELSE IF (P. GT. PG) THEN
                  IREG3=2
               ELSE 
                  IREG3=3
               ENDIF
            ELSE IF (P .GT. PG) THEN
               IREG3=5
            ELSE
               IREG3=6
            ENDIF
         ENDIF
      ELSE
         IREG3=7
      ENDIF
C
      END            
C
C***********************************************************************
        DOUBLEPRECISION FUNCTION VEST3(P,T,IREG3)
C***********************************************************************
        implicit double precision(a-h,o-z)
C
      COMMON/CSUB2/TC,PC,DC
C
        parameter(vc=1.d0/322.d0)
        parameter(dtdp=3.727888d0)
c parameters for region 3a
        parameter(a01= 0.748874448725080D-02)
        parameter(a02= 0.348351238844209D-03)
        parameter(a03=-0.118427102787648D-03)
        parameter(a04=-0.332156380721166D-01)
        parameter(a05=-0.797144480049465D-03)
        parameter(a06= 0.608359392259313D-01)
        parameter(a07= 0.617493752916276D-03)
        parameter(a08=-0.564778221803544D-01)
        parameter(a09=-0.161628813019544D-03)
        parameter(a10=-0.876842774061156D-02)
        parameter(a11= 0.263271308078056D-01)
        parameter(a12= 0.172887471616688D-01)
        parameter(a13=-0.490643103428227D-02)
        parameter(a14=-0.106812857513379D-01)
        parameter(a15= 0.153431224177324D-02)
        parameter(a16=-0.383218936572771D-08)
        parameter(a17=-0.455907876060086D-04)
        parameter(a18=-0.251644090693395D-11)
c parameters for region 3b
        parameter(b01=-0.487687959756292D-01)
        parameter(b02= 0.135336193587336D-02)
        parameter(b03= 0.413013399556249d0  )
        parameter(b04= 0.122945704230431d0  )
        parameter(b05=-0.141140213612559D-03)
        parameter(b06= 0.178357790645547D-01)
        parameter(b07= 0.395174104863528D-07)
        parameter(b08=-0.585398228522495d0  )
        parameter(b09= 0.594938907295817D-04)
        parameter(b10=-0.119057019174713d0  )
        parameter(b11=-0.194570794816001D-02)
        parameter(b12= 0.427529172119353D-01)
        parameter(b13=-0.760561375637742D-01)
        parameter(b14= 0.299546836134716d0  )
        parameter(b15= 0.139391657221250D-01)
        parameter(b16=-0.278468126229603D-01)
        parameter(b17= 0.195041855369227D-02)
        parameter(b18= 0.375831122504734D-04)
        parameter(b19= 0.125685384163741D-03)
        parameter(b20=-0.762583183784896D-02)
        parameter(b21=-0.938937201479048D-06)
        parameter(b22= 0.134876477124427D-01)
        parameter(b23= 0.634557820845727D-02)
        parameter(b24=-0.291792573898258D-01)
        parameter(b25=-0.136017045350123D-05)
        parameter(b26= 0.454845318075243D-07)
        parameter(b27=-0.239512088942747D-09)
c parameters for region 3c
        parameter(c01=-0.151226672652382D-01)
        parameter(c02= 0.219075237766159D-03)
        parameter(c03= 0.114144924756274D-04)
        parameter(c04=-0.470815410341398D-02)
        parameter(c05= 0.105510662596481D-03)
        parameter(c06= 0.487932009131791D0  )
        parameter(c07= 0.424198281757227D-03)
        parameter(c08= 0.305426466180436D-01)
        parameter(c09=-0.174690467895005D-04)
        parameter(c10= 0.509486478795057D-06)
        parameter(c11= -1.43708991982910D0  )
        parameter(c12=-0.526750160303121D-04)
        parameter(c13=  1.22116000890535D0  )
        parameter(c14= 0.104163340234817D-06)
        parameter(c15= 0.218460997189951D-02)
        parameter(c16=-0.152017319222412D-01)
        parameter(c17=-0.294106550573793D-01)
        parameter(c18= 0.165284534427183D-01)
        parameter(c19=-0.335234582911578D0  )
        parameter(c20=-0.212663865441498D-02)
        parameter(c21= 0.753543651141502D-01)
        parameter(c22=-0.769927079971342D-01)
        parameter(c23=-0.247039860992736D-01)
        parameter(c24=-0.223745526548978D-02)
        parameter(c25=-0.113782324304171D-01)
        parameter(c26= 0.122241653100711D-01)
        parameter(c27=-0.163238458626065D-06)
        parameter(c28= 0.580172883857322D-09)
        parameter(c29=-0.419492118766324D-11)
        parameter(c30=-0.660478916724586D-06)
        parameter(c31=-0.119419038095450D-05)
        parameter(c32= 0.310771050242050D-01)
c parameters for region 3d
        parameter(d01= 0.664933791260178D-01)
        parameter(d02= 0.571985898251817D-06)
        parameter(d03=-0.274711031859796D0  )
        parameter(d04=-0.252680942664202D0  )
        parameter(d05=-0.392289014485859D-06)
        parameter(d06= 0.256374705633001D-10)
        parameter(d07=  1.04023037508896D0  )
        parameter(d08= 0.643907516287541D-03)
        parameter(d09= 0.389907326497054D0  )
        parameter(d10= -1.60832918484839D0  )
        parameter(d11= 0.713823865906627D-03)
        parameter(d12=-0.303433314767063D0  )
        parameter(d13=  1.25833184484538D0  )
        parameter(d14=-0.792816174523651D-04)
        parameter(d15= 0.118639803204261D0  )
        parameter(d16=-0.496046407723234D0  )
        parameter(d17=-0.185997001981635D-01)
        parameter(d18= 0.785791469799311D-01)
c parameters for region 3e
        parameter(e01= 0.301448430808593D-01)
        parameter(e02=-0.335916193888413D-03)
        parameter(e03=-0.112702460700323D0  )
        parameter(e04=-0.877499173794533D-20)
        parameter(e05= 0.187152530917710D-03)
        parameter(e06= 0.356911241995526D-02)
        parameter(e07=-0.418986486825425D-11)
        parameter(e08= 0.322726674818840D-06)
        parameter(e09=-0.120527379090314D-10)
        parameter(e10=-0.817132957960820D-13)
        parameter(e11=-0.289296115834392D-07)
        parameter(e12=-0.194260348835721D-06)
        parameter(e13=-0.813610292866497D-09)
        parameter(e14= 0.391792707971363D-13)
        parameter(e15=-0.519162117274822D-07)
        parameter(e16=-0.538255397523665D-06)
        parameter(e17= 0.200848495263495D-07)
        parameter(e18=-0.721941889977446D-02)
        parameter(e19=-0.313609300142694D-02)
        parameter(e20= 0.260485012629641D-04)
        parameter(e21=-0.370031971083042D-03)
        parameter(e22= 0.759164825488741D-07)
        parameter(e23=-0.169209023050985D-09)
        parameter(e24= 0.396910783770869D-01)
        parameter(e25= 0.127931680641201D-01)
        parameter(e26= 0.191456960283807D-06)
        parameter(e27= 0.766425815924100D-03)
        parameter(e28=-0.431282885175170D-05)
        parameter(e29= 0.225545818343096D-04)
        parameter(e30= 0.488412051812406D-07)
        parameter(e31=-0.356918730340587D-09)
        parameter(e32= 0.303625591962761D-03)
        parameter(e33= 0.124536067659580D-12)
        parameter(e34=-0.550701260682904D-02)
        parameter(e35= 0.252719052207758D-11)
        parameter(e36= 0.133491926757520D-12)
        parameter(e37=-0.607043036830454D-02)
        parameter(e38= 0.287671539256386D-05)
        parameter(e39=-0.958920650759504D-04)
        parameter(e40= 0.102650393925122D-07)
        parameter(e41=-0.155995837253683D-15)
c parameters for region 3f
        parameter(f01= 0.329890852403526D-04)
        parameter(f02= 0.379781955709120D-04)
        parameter(f03= 0.252360667718127D-06)
        parameter(f04=-0.407113420864393D-03)
        parameter(f05=-0.357933830043717D-04)
        parameter(f06=-0.245865969502991D-03)
        parameter(f07= 0.525013249033141D-09)
        parameter(f08= 0.321621088787326D-02)
        parameter(f09= 0.439053317238974D-03)
        parameter(f10=-0.587768799763869D-07)
        parameter(f11= 0.262763491213355D-05)
        parameter(f12=-0.214747822942469D-02)
        parameter(f13=-0.635909964040088D-04)
        parameter(f14=-0.115529419916092D-06)
        parameter(f15= 0.785448130272979D-03)
        parameter(f16=-0.259174703895765D-15)
        parameter(f17= 0.184266268449228D-10)
c parameters for region 3g
        parameter(g01= 0.125537070183712D-02)
        parameter(g02= 0.130245527248376D-02)
        parameter(g03= 0.103367194207180D-01)
        parameter(g04=-0.254795720214314D-01)
        parameter(g05=-0.185955960512067D-03)
        parameter(g06=-0.960082399513164D-02)
        parameter(g07= 0.474944869074855D0  )
        parameter(g08=-0.503420527214133D-02)
        parameter(g09=-0.687909934564732D0  )
        parameter(g10=  30.7310678953686D0  )
        parameter(g11=-0.168658389645091D-03)
        parameter(g12=  654.223452156635D0  )
        parameter(g13= -320317.604761443D0  )
        parameter(g14=  197.974246206705D0  )
        parameter(g15=  11416042249924.4D0  )
        parameter(g16= -9156698413.12590D0  )
        parameter(g17=  2783392577409.60D0  )
        parameter(g18= -281900719117892.D0  )
        parameter(g19= 0.119502836257688D+17)
        parameter(g20=-0.183657231751509D+18)
C
      GOTO (100,200,300,400,500,600,700) IREG3            
c region 3a (t<tc, p<24 MPa, liquid)
  100   pr=p-pc
        tr=t-tc
        z8=(dtdp+.05d0*pr)*pr-(1.d0-.05d0*tr)*tr
c to avoid negative roots outside of range of validity
        if (z8.lt.0.d0) then
          z4=-1.d0*sqrt(-1.d0*z8)
          z2=-1.d0*sqrt(-1.d0*z4)
          z1=-1.d0*sqrt(-1.d0*z2)
        else
          z4=sqrt(z8)
          z2=sqrt(z4)
          z1=sqrt(z2)
        endif
        vpt3n=a03*z1+z4*(a10+z1*(a12+z1*(a14+z2*(a15+z4*a17))))
     +  +tr*(a01+z1*(a04+z1*(a06+z1*(a08+z1*(a11+z1*a13))))
     +  +tr*(a02+z1*(a05+z1*(a07+z1*(a09+tr*z1*z4*(a16+tr*z4*a18))))))
        vest3=vpt3n+vc
C
      GOTO 999
C      
c region 3b (t>tc, v<~vc, p<24 MPa)
  200 pr=p-pc
        p2=pr*pr
        p4=p2*p2
        tr=t-tc
        t2=tr*tr
        z8=(dtdp+.02d0*pr)*pr-(1.d0-.02d0*tr)*tr
c to avoid negative roots outside of range of validity
        if (z8.lt.0.d0) then
          z4=-1.d0*sqrt(-1.d0*z8)
          z2=-1.d0*sqrt(-1.d0*z4)
          z1=-1.d0*sqrt(-1.d0*z2)
        else
          z4=sqrt(z8)
          z2=sqrt(z4)
          z1=sqrt(z2)
        endif
        vpt3n=p2*b03+t2*(b01+pr*b02)
     +  +z1*(pr*b06+p2*b08+tr*(p4*b09+tr*(b04+tr*(b05+t2*pr*b07)))
     +  +z1*(t2*b10
     +  +z1*(b11+pr*b13+p2*b14+t2*b12
     +  +z1*(tr*b15
     +  +z1*(tr*pr*b16
     +  +z1*(b17
     +  +z2*(p2*(pr*b19+t2*b18)
     +  +z2*(tr*b20+t2*p4*pr*b21
     +  +z2*(b22+pr*b24+tr*(b23+t2*p2*(b25+p4*(b26+t2*b27))) )))))))))
        vest3=vpt3n+vc
C       
      GOTO 999
C       
c region 3c (v>~vc, pc<p<24 MPa)
  300 pr=p-pc
        p2=pr*pr
        tr=t-tc
        t2=tr*tr
        t3=t2*tr
        t4=t2*t2
        z8=(.2d0*pr-dtdp)*pr+(1.d0+.02d0*tr)*tr
c to avoid negative roots outside of range of validity
        if (z8.lt.0.d0) then
          z4=-1.d0*sqrt(-1.d0*z8)
          z2=-1.d0*sqrt(-1.d0*z4)
          z1=-1.d0*sqrt(-1.d0*z2)
        else
          z4=sqrt(z8)
          z2=sqrt(z4)
          z1=sqrt(z2)
        endif
        vpt3n=t2*c01+t3*c02+t4*c03+pr*(c04+t3*c05+pr*(c06+p2*c07))
     +  +z1*(t2*c08+t4*(c09+tr*c10)+p2*c11
     +  +z1*(pr*t3*c12+p2*(c13+t4*c14)
     +  +z1*(c15+t2*c16+pr*(c17+tr*c18)+p2*(c19+pr*c20)
     +  +z2*(tr*c21
     +  +z1*(tr*c22+pr*c23
     +  +z2*(p2*c24
     +  +z2*(c25
     +  +tr*(c26+t3*c27
     +  +z2*(t4*(c28+tr*c29)+p2*(tr*c30+pr*c31) ))) ))) )))
        vest3=vpt3n+vc
C       
        GOTO 999
C       
c region 3d (t>tsat, p<pc)
  400 pr=p-pc
        tr=t-tc
        t2=tr*tr
        z8=(.5d0*pr-dtdp)*pr+(1.d0+.01d0*tr)*tr
c to avoid negative roots outside of range of validity
        if (z8.lt.0.d0) then
          z4=-1.d0*sqrt(-1.d0*z8)
          z2=-1.d0*sqrt(-1.d0*z4)
          z1=-1.d0*sqrt(-1.d0*z2)
        else
          z4=sqrt(z8)
          z2=sqrt(z4)
          z1=sqrt(z2)
        endif
        vpt3n=z2*(d08+z1*(d11+z1*d14))
     +  +tr*(d01+z1*(d04+z1*(d09+z1*(d12+z1*(d15+z1*d17))))
     +  +t2*(d02+z1*(d05+t2*d06)))
     +  +pr*(d03+z1*(d07+z1*(d10+z1*(d13+z1*(d16+z1*d18)))))
        vest3=vpt3n+vc
C       
        GOTO 999
C       
c region 3e (v<~vc, 24 MPa<p<40 MPa)
  500 pr=p-pc
        p2=pr*pr
        p4=p2*p2
        p5=p4*pr
        tr=t-tc
        t2=tr*tr
        t3=t2*tr
        z8=(dtdp+.02d0*pr)*pr+(.005d0*tr-1.d0)*tr
c to avoid negative roots outside of range of validity
        if (z8.lt.0.d0) then
          z4=-1.d0*sqrt(-1.d0*z8)
          z2=-1.d0*sqrt(-1.d0*z4)
          z1=-1.d0*sqrt(-1.d0*z2)
        else
          z4=sqrt(z8)
          z2=sqrt(z4)
          z1=sqrt(z2)
        endif
        vpt3n=pr*e03+tr*e01+t2*(e02+t2*t2*p5*e04)
     +  +z1*(t2*e05
     +  +z1*(pr*e06+t2*p5*e07
     +  +z1*(tr*p4*e11+t3*(e08+tr*p2*pr*e10+t2*e09)
     +  +z1*(p5*e15+t3*(e12+tr*e13+t3*e14)
     +  +z1*(p4*(e16+tr*e17)
     +  +z1*(e18+tr*(e19+pr*e21+p5*e23+tr*(e20+p2*e22))
     +  +z2*(e24+pr*(e25+pr*e27+p2*e29+tr*(pr*e28+p2*e30+tr*(e26+p2*e31
     +          )))
     +  +z2*(pr*e34+tr*(e32+p5*pr*e36+tr*(p4*e35+t3*e33))
     +  +z2*(e37+p2*e39+p4*e40+tr*pr*(e38+tr*p5*e41) ))) ))) )))
        vest3=vpt3n+vc
C       
        GOTO 999
C       
c region 3f (v>~vc, 24 MPa<p<40 MPa)
  600   pr=p-pc
        p2=pr*pr
        tr=t-tc
        z4=(.1d0*pr-dtdp)*pr+(1.d0+.002d0*tr)*tr
c to avoid negative roots outside of range of validity
        if (z4.lt.0.d0) then
          z2=-1.d0*sqrt(-1.d0*z4)
          z1=-1.d0*sqrt(-1.d0*z2)
        else
          z2=sqrt(z4)
          z1=sqrt(z2)
        endif
        z3=z1*z2
        vpt3n=z3*(f08+z2*(f12+p2*f14+z1*f15))
     +  +tr*(f01+z1*(f04+z1*(f06+z1*(f09+z2*(f13+z1*p2*pr*f17))))
     +  +tr*(f02+z1*(f05+z3*f11)
     +  +tr*(f03+z3*f10
     +  +tr*z2*(f07+z4*p2*f16))))
        vest3=vpt3n+vc
C       
        GOTO 999
C       
c region 3g (40 MPa<p<100 MPa)
  700   pr=p*(1.d0/103.d0)-1.d0
        p2=pr*pr
        p4=p2*p2
        tr=t*(1.d0/600.d0)-1.d0
        t2=tr*tr
        t3=t2*tr
        t4=t2*t2
        t6=t3*t3
        p2t2=p2*t2
        Z8=c32*p2
        vpt3n=g01+tr*(g02+t2*(g03+t3*g04))
     +  +pr*(g05+t2*g06
     +  +pr*(tr*t4*g07
     +  +pr*(tr*(g08+t3*(g09+t4*g10))
     +  +p2*(g11
     +  +p2*pr*(t6*g12
     +  +p2*t4*(t6*g13
     +  +p4*(g14
     +  +p2*t6*(t6*t6*g15
     +  +p4*p4*(g16
     +  +p2t2*(g17
     +  +p2t2*(g18
     +  +p2t2*(g19+p2t2*g20))) ))) ))) )))
      vest3=vpt3n
C     
      GOTO 999
C     
  999 CONTINUE
C
        END
C
C************************************************************************
      DOUBLEPRECISION FUNCTION DITER3(P,T,DEST,EPS)
C************************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
      double precision NULLP3N
      EXTERNAL NULLP3N
c
      COMMON/CSUB2/TC,PC,DC
c
      d1v=dest*0.999d0
        d2v=dest*1.001d0
        If(dest .le. dc) then
         d1=d1v
           if (d2v .le. dc) then 
              d2=d2v
           else
              d2=dc
           endif
        else
           d2=d2v
           if (d1v .ge. dc) then 
              d1=d1v
           else
              d1=dc
           endif
        endif
c
      call WNPT3(d1,d2,nullp3n,p,t,eps,x,ix)
c
      if( ix .le. 0 ) then
        DITER3=x
      else
        DITER3=0.d0
      end if
      end
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NULLP3N(D,T,P)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
c
      V=1.D0/D
      NULLP3N=PVT3N(V,T)-P
c
      END
C
C*******************************************************************
      SUBROUTINE WNPT3(XA,XB,F,P,T,EPS,X,IX)
C*******************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      EXTERNAL F
      X1=XA
      F1=F(X1,T,P)
      X3=XB
      F3=F(X3,T,P)
      IX=0
C++ SCHLEIFE ZUR ERMITTLUNG DER NULLSTELLE +++++++++++++++++++++++++++
      DO 100 I=1,40
      IF( F1 .NE. F3 ) THEN
       X=X1+(X3-X1)*F1/(F1-F3)
      ELSE
       IX=3
       GOTO 999
      END IF
      IF(X .LT. 0.D0) X=(X1+X3)/2.D0
      IF(DABS(X) .LT. 1.D-8) THEN
        IF(DABS(X-X1).LT.EPS) RETURN
      ELSE
        IF(DABS((X-X1)/X).LT.EPS) RETURN
      END IF
      F2=F(X,T,P)
      X2=X1-(X1-X3)/2.D0
      P1=F2*F1
      P3=F2*F3
      IF ((P1 .LT. 0.D0) .AND. (P3 .LT. 0.D0)) THEN
        P1=DABS(P1)
        P3=DABS(P3)
      ENDIF
      IF(P1 .LE. P3) THEN
        X3=X1
        F3=F1
      END IF
      X1=X
      F1=F2
      IF((X2-X3)*(X2-X1) .GE. 0.D0) GOTO 100
      IF ((DABS(F1/F3)) .LE. (0.3D0)) THEN
         X=X1+DABS(F1/F3)*1.5D0*(X3-X1)
         IREM=1
      ELSE IF ((DABS(F1/F3)) .GE. (3.D0)) THEN
         IREM=2
      ELSE
         X=(X1+X3)/2.D0
         IREM=0
      ENDIF
      F2=F(X,T,P)
      IF ((((F2*F1) .GE. 0.D0) .AND. (IREM .EQ. 1)) .OR.
     1    (((F2*F3) .GE. 0.D0) .AND. (IREM .EQ. 2))) THEN
         X=(X1+X3)/2.D0
         F2=F(X,T,P)
      ENDIF
      IF(F2*F1 .LE. F2*F3) THEN
        X3=X1
        F3=F1
      END IF
      X1=X
      F1=F2
 100  CONTINUE
C++ ENDE DER SCHLEIFE ++++++++++++++++++++++++++++++++++++++++++++++++
      IX=1
 999  RETURN
      END
C
C***********************************************************************
      SUBROUTINE REGSPH(P,H,IREG)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      COMMON/CSUB2/TC,PC,DC
      COMMON/CGR/PGR,TGR13
C
      double precision HOLD, POLD
      SAVE HOLD, POLD
        INTEGER IREGOLD
        SAVE IREGOLD
      data HOLD, POLD, IREGOLD / -1.D-6, -1.D-6, -1 /
C
      IF ((DABS(H-HOLD).LT. 1.D-6) .AND.
     *   (DABS(P-POLD).LT. 1.D-6)) THEN
         IREG=IREGOLD
         GOTO 999
      END IF
C
      IF ((P .GT. 100.D0) .OR. (P .LT. 5.D-4) .OR.
     1    (H .LT. -0.05D0) .OR. (H .GT. 7.4D3)) THEN
         IREG=0
         GOTO 999
      ENDIF
C
        IREG=0
C
      IF (P .LE. PGR) THEN
         TG=TSATPN(P)
         H1G=HPT1N(P,TG)
         IF (H .LE. H1G) THEN
            TN=TPH1N(P,H)
            IF ((TN .GE. 273.125D0) .AND. (TN .LT. 273.15D0)) THEN
               TN=273.15D0
            ENDIF
            IF (TN .GE. 273.15D0) THEN
               IREG=1
            ELSE
               IREG=0
            ENDIF
         ELSE
            H2G=HPT2N(P,TG)
            IF (H .LT. H2G) THEN
               IREG=9
            ELSE IF (H .LE. 4080.D0) THEN
               IREG=2
            ELSE
               HG2=HPT2N(P,1073.1501D0)
                 IF (H .LE. HG2) THEN
                    IREG=2
               ELSE IF (P .GT. 10.D0) THEN
                  IREG=0
               ELSE IF (H .LT. 7374.751711911952D0) THEN
                  IREG=5
               ELSE IF (H .LT. 7376.980306192378D0) THEN
                  TB=TPH5N(P,H)
                  IF (TB .LE. 2273.15D0) THEN
                     IREG=5
                  ELSE
                     IREG=0
                  ENDIF
               ELSE 
                  IREG=0
               ENDIF
            ENDIF
         ENDIF
      ELSE IF (H .GE. 2563.592003888252D0) THEN
         IF (H .GT. 4090D0) THEN
            IREG=0
         ELSE
            TG=FB23P(P)
              HG23=HPT2N(P,TG)
            IF (H .LT. HG23) THEN
                 IREG=3
              ELSE
                 HG20=HPT2N(P,1073.175D0)
                 IF (H .LE. HG20) THEN
                    IREG=2
                 ELSE
                    IREG=0
                 ENDIF
              ENDIF
           ENDIF
      ELSE
         HG13=HPT1N(P,623.15D0)
           IF (H .GT. HG13) THEN
              IREG=3
           ELSE
              HG01=HPT1N(P,273.125D0)
              IF (H .GE. HG01) THEN
                 IREG=1
              ELSE
                 IREG=0
              ENDIF
           ENDIF
        ENDIF
      IF ((IREG .EQ. 3) .AND. (P .GT. PGR) .AND. (P .LT. PC)) THEN
         TS=TSATPN(P)
         CALL FSATP(DV,DL,TOUT,TS,P)
         V1=1.D0/DL
         V2=1.D0/DV
         H1=HVT3N(V1,TS)
         H2=HVT3N(V2,TS)
         IF ((H .GT. H1) .AND. (H .LT. H2)) THEN
            IREG=9
         ENDIF
      ENDIF
C
999   CONTINUE
C      
      HOLD=H
      POLD=P
      IREGOLD=IREG
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION TPH(P,H,IREG)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      IF (IREG .EQ. 0) THEN
         CALL REGSPH(P,H,IREG)
      ENDIF
C  
      IF (IREG .EQ. 1) THEN
         TPH=TPH1N(P,H)
         IF ((TPH .GE. 273.125D0) .AND. (TPH .LT. 273.15D0)) THEN
            TPH=273.15D0
         ENDIF
      ELSE IF (IREG .EQ. 2) THEN
         TPH=TPH2N(P,H)
         IF ((TPH .GT. 1073.15D0) .AND. (TPH .LE. 1073.175D0)) THEN
            TPH=1073.15D0
         ENDIF
      ELSE IF (IREG .EQ. 3) THEN
         TPH=TPH3N(P,H)
      ELSE IF (IREG .EQ. 5) THEN
         TPH=TPH5N(P,H)
      ELSE IF (IREG .EQ. 9) THEN
         TPH=TSATPN(P)
      ELSE
         TPH=0.D0
      ENDIF
C
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION TPH3N(P,H)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      EPS=1.D-9
C
      TE=TPH3V(P,H)
C
      TPH3N=TITH3(P,H,TE,EPS)
C
      END      
C
C************************************************************************
      DOUBLE PRECISION FUNCTION TPH3V(P,H)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      COMMON/CSUB2/TC,PC,DC
C
        parameter(hc=2087.54684511650d0)
        parameter(a1=  646.836937239618D0  )
        parameter(a2=  62.7217723183194D0  )
        parameter(a3= -18.0176207184787D0  )
        parameter(a4=  2.01708278241311D0  )
        parameter(a5=  247.531572315595D0  )
        parameter(a6= -41.1717709139789D0  )
        parameter(a7= -150.208861198343D0  )
        parameter(a8=  1431.77901876243D0  )
        parameter(a9= -862.734816079672D0  )
        parameter(a0=  49.9470109089612D0  )
        hr=h*(1.d0/hc)-1.d0
        pr=p*(1.d0/pc)-1.d0
        TPH3V=a1+pr*(a2+pr*(a3+pr*a4))
     +  +hr*(pr*(a5+pr*a6)
     +  +hr*(a7+hr*(a8+pr*(a9+pr*pr*a0))))
C
        END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION TITH3(P,H,TE,EPS)
C************************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
      double precision NULLH3N
      EXTERNAL NULLH3N
C
      IF (TE .LT. 629.D0) THEN
         T1=623.15D0
      ELSE
         T1=TE*0.996D0
      ENDIF
C
      TG=FB23P(P)
      IF (TE .GT. (TG/1.01D0)) THEN
         T2=TG
      ELSE   
         T2=TE*1.004D0
      ENDIF
C
      CALL WNPH3(T1,T2,NULLH3N,H,P,EPS,X,IX)
C
      IF( IX .LE. 0 ) THEN
        TITH3=X
      ELSE
        TITH3=0.d0
      END IF
C      
      END
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NULLH3N(T,P,H)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      VE=VPT3N(P,T)
      NULLH3N=HVT3N(VE,T)-H
C
      END
C
C*******************************************************************
      SUBROUTINE WNPH3(XA,XB,F,P,T,EPS,X,IX)
C*******************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      EXTERNAL F
      X1=XA
      F1=F(X1,T,P)
      X3=XB
      F3=F(X3,T,P)
      IX=0
C++ SCHLEIFE ZUR ERMITTLUNG DER NULLSTELLE +++++++++++++++++++++++++++
      DO 100 I=1,40
      IF( F1 .NE. F3 ) THEN
       X=X1+(X3-X1)*F1/(F1-F3)
      ELSE
       IX=3
       GOTO 999
      END IF
      IF(X .LT. 0.D0) X=(X1+X3)/2.D0
      IF(DABS(X) .LT. 1.D-8) THEN
        IF(DABS(X-X1).LT.EPS) RETURN
      ELSE
        IF(DABS((X-X1)/X).LT.EPS) RETURN
      END IF
      F2=F(X,T,P)
      X2=X1-(X1-X3)/2.D0
      P1=F2*F1
      P3=F2*F3
      IF ((P1 .LT. 0.D0) .AND. (P3 .LT. 0.D0)) THEN
        P1=DABS(P1)
        P3=DABS(P3)
      ENDIF
      IF(P1 .LE. P3) THEN
        X3=X1
        F3=F1
      END IF
      X1=X
      F1=F2
      IF((X2-X3)*(X2-X1) .GE. 0.D0) GOTO 100
c      X=(X1+X3)/2.D0
      IF ((DABS(F1/F3)) .LE. (0.3D0)) THEN
         X=X1+DABS(F1/F3)*1.5D0*(X3-X1)
         IREM=1
      ELSE IF ((DABS(F1/F3)) .GE. (3.D0)) THEN
         IREM=2
      ELSE
         X=(X1+X3)/2.D0
         IREM=0
      ENDIF
      F2=F(X,T,P)
      IF ((((F2*F1) .GE. 0.D0) .AND. (IREM .EQ. 1)) .OR.
     1    (((F2*F3) .GE. 0.D0) .AND. (IREM .EQ. 2))) THEN
         X=(X1+X3)/2.D0
         F2=F(X,T,P)
      ENDIF
      IF(F2*F1 .LE. F2*F3) THEN
        X3=X1
        F3=F1
      END IF
      X1=X
      F1=F2
 100  CONTINUE
C++ ENDE DER SCHLEIFE ++++++++++++++++++++++++++++++++++++++++++++++++
      IX=1
 999  RETURN
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION TPH5N(P,H)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision NULLH5N
      EXTERNAL NULLH5N
C
      EPS=1.D-6
C
      IF (H .LT. 5000.D0) THEN
         T1=1073.15D0
         T2=1573.15D0
      ELSE IF (H .LT. 6500.D0) THEN
         T1=1323.15D0
         T2=2073.15D0
      ELSE
         T1=1823.15D0
         T2=2273.15D0
      ENDIF
C
      CALL WNPH3(T1,T2,NULLH5N,H,P,EPS,X,IX)
C      
      TPH5N=X
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NULLH5N(T,P,H)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      NULLH5N=HPT5N(P,T)-H
C
      END
C
C***********************************************************************
      SUBROUTINE REGSPS(P,S,IREG)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      COMMON/CSUB2/TC,PC,DC
      COMMON/CGR/PGR,TGR13
C
      double precision POLD, SOLD
      SAVE POLD, SOLD
        INTEGER IREGOLD
        SAVE IREGOLD
      data POLD, SOLD, IREGOLD / -1.D+0, -1.D-5, -1 /
C
      IF ((DABS(S-SOLD).LT. 1.D-6) .AND.
     *   (DABS(P-POLD).LT. 1.D-6)) THEN
         IREG=IREGOLD
         GOTO 999
      END IF
C
      IF ((P .GT. 100.D0) .OR. (P .LT. 5.D-4) .OR.
     1    (S .LT. -0.5D0).OR. (S .GT. 14.D0)) THEN
         IREG=0
         GOTO 999
      ENDIF
C
        IREG=0
C
      IF (P .LE. PGR) THEN
         TG=TSATPN(P)
         S1G=SPT1N(P,TG)
         IF (S .LE. S1G) THEN
            TN=TPS1N(P,S)
            IF ((TN .GE. 273.125D0) .AND. (TN .LT. 273.15D0)) THEN
               TN=273.15D0
            ENDIF
            IF (TN .GE. 273.15D0) THEN
               IREG=1
            ELSE
               IREG=0
            ENDIF
         ELSE
            S2G=SPT2N(P,TG)
            IF (S .LT. S2G) THEN
               IREG=9
            ELSE IF (S .LE. 7.15D0) THEN
               IREG=2
            ELSE
               SG2=SPT2N(P,1073.1501D0)
               IF (S .LE. SG2) THEN
                  IREG=2
               ELSE IF (P .GT. 10.D0) THEN
                  IREG=0
               ELSE IF (S .LT. 9.423909357921552D0) THEN
                  IREG=5
               ELSE IF (S .LT. 13.99764756652987D0) THEN
                  TB=TPS5N(P,S)
                  IF (TB .LE. 2273.2D0) THEN
                     IREG=5
                  ELSE
                     IREG=0
                  ENDIF
               ELSE 
                  IREG=0
               ENDIF
            ENDIF
         ENDIF
      ELSE IF (S .GE. 3.778281339544279D0) THEN
         IF (S .GT. 7.2D0) THEN
            IREG=0
         ELSE
            TG=FB23P(P)
              SG23=SPT2N(P,TG)
              IF (S .LT. SG23) THEN
                 IREG=3
              ELSE
                 SG20=SPT2N(P,1073.175D0)
                 IF (S .LE. SG20) THEN
                    IREG=2
                 ELSE
                    IREG=0
                 ENDIF
              ENDIF
         ENDIF
      ELSE
         SG13=SPT1N(P,623.15D0)
           IF (S .GT. SG13) THEN
              IREG=3
           ELSE
              SG01=SPT1N(P,273.125D0)
              IF (S .GE. SG01) THEN
                 IREG=1
              ELSE
                 IREG=0
              ENDIF
           ENDIF
      ENDIF      
      IF ((IREG .EQ. 3) .AND. (P .GT. PGR) .AND. (P .LT. PC)) THEN
         TS=TSATPN(P)
         CALL FSATP(DV,DL,TOUT,TS,P)
         V1=1.D0/DL
         V2=1.D0/DV
         S1=SVT3N(V1,TS)
         S2=SVT3N(V2,TS)
         IF ((S .GT. S1) .AND. (S .LT. S2)) THEN
            IREG=9
         ENDIF
      ENDIF
C
999   CONTINUE
C      
      POLD=P
      SOLD=S
      IREGOLD=IREG
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION TPS(P,S,IREG)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      IF (IREG .EQ. 0) THEN
         CALL REGSPS(P,S,IREG)
      ENDIF
C  
      IF (IREG .EQ. 1) THEN
         TPS=TPS1N(P,S)
         IF ((TPS .GE. 273.125D0) .AND. (TPS .LT. 273.15D0)) THEN
            TPS=273.15D0
         ENDIF
      ELSE IF (IREG .EQ. 2) THEN
         TPS=TPS2N(P,S)
         IF ((TPS .GT. 1073.15D0) .AND. (TPS .LE. 1073.175D0)) THEN
            TPS=1073.15D0
         ENDIF
      ELSE IF (IREG .EQ. 3) THEN
         TPS=TPS3N(P,S)
      ELSE IF (IREG .EQ. 5) THEN
         TPS=TPS5N(P,S)
      ELSE IF (IREG .EQ. 9) THEN
         TPS=TSATPN(P)
      ELSE
         TPS=0.D0
      ENDIF
C
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION TPS3N(P,S)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      EPS=1.D-9
C
      TE=TPS3V(P,S)
C
      TPS3N=TITS3(P,S,TE,EPS)
C
      END      
C
C************************************************************************
      DOUBLE PRECISION FUNCTION TPS3V(P,S)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      COMMON/CSUB2/TC,PC,DC
C
        parameter(sc=4.41202148223476d0)
        parameter(a1=  646.350957423454D0  )
        parameter(a2=  67.2054900478099D0  )
        parameter(a3= -15.4539626100284D0  )
        parameter(a4=  1.64063096179835D0  )
        parameter(a5=  339.504918981144D0  )
        parameter(a6= -80.8665059602751D0  )
        parameter(a7=  8.49001610510301D0  )
        parameter(a8= -260.117974798217D0  )
        parameter(a9=  3252.66519415985D0  )
        parameter(a0=  4686.59047913868D0  )
        sr=s*(1.d0/sc)-1.d0
        pr=p*(1.d0/pc)-1.d0
        TPS3V=a1+pr*(a2+pr*(a3+pr*a4))
     +  +sr*(pr*(a5+pr*(a6+pr*a7))
     +  +sr*(a8+sr*(a9+sr*pr*a0)))
C       
        END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION TITS3(P,S,TE,EPS)
C************************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
      double precision NULLS3N
      EXTERNAL NULLS3N
C
      IF (TE .LT. 629.D0) THEN
         T1=623.15D0
      ELSE
         T1=TE*0.995D0
      ENDIF
C
      TG=FB23P(P)
      IF (TE .GT. (TG/1.01D0)) THEN
         T2=TG
      ELSE   
         T2=TE*1.005D0
      ENDIF
C
      CALL WNPS3(T1,T2,NULLS3N,S,P,EPS,X,IX)
C
      IF( IX .LE. 0 ) THEN
        TITS3=X
      ELSE
        TITS3=0.d0
      END IF
C      
      END
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NULLS3N(T,P,S)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      VE=VPT3N(P,T)
      NULLS3N=SVT3N(VE,T)-S
C
      END
C
C*******************************************************************
      SUBROUTINE WNPS3(XA,XB,F,P,T,EPS,X,IX)
C*******************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      EXTERNAL F
      X1=XA
      F1=F(X1,T,P)
      X3=XB
      F3=F(X3,T,P)
      IX=0
C++ SCHLEIFE ZUR ERMITTLUNG DER NULLSTELLE +++++++++++++++++++++++++++
      DO 100 I=1,40
      IF( F1 .NE. F3 ) THEN
       X=X1+(X3-X1)*F1/(F1-F3)
      ELSE
       IX=3
       GOTO 999
      END IF
      IF(X .LT. 0.D0) X=(X1+X3)/2.D0
      IF(DABS(X) .LT. 1.D-8) THEN
        IF(DABS(X-X1).LT.EPS) RETURN
      ELSE
        IF(DABS((X-X1)/X).LT.EPS) RETURN
      END IF
      F2=F(X,T,P)
      X2=X1-(X1-X3)/2.D0
      P1=F2*F1
      P3=F2*F3
      IF ((P1 .LT. 0.D0) .AND. (P3 .LT. 0.D0)) THEN
        P1=DABS(P1)
        P3=DABS(P3)
      ENDIF
      IF(P1 .LE. P3) THEN
        X3=X1
        F3=F1
      END IF
      X1=X
      F1=F2
      IF((X2-X3)*(X2-X1) .GE. 0.D0) GOTO 100
c      X=(X1+X3)/2.D0
      IF ((DABS(F1/F3)) .LE. (0.3D0)) THEN
         X=X1+DABS(F1/F3)*1.5D0*(X3-X1)
         IREM=1
      ELSE IF ((DABS(F1/F3)) .GE. (3.D0)) THEN
         IREM=2
      ELSE
         X=(X1+X3)/2.D0
         IREM=0
      ENDIF
      F2=F(X,T,P)
      IF ((((F2*F1) .GE. 0.D0) .AND. (IREM .EQ. 1)) .OR.
     1    (((F2*F3) .GE. 0.D0) .AND. (IREM .EQ. 2))) THEN
         X=(X1+X3)/2.D0
         F2=F(X,T,P)
      ENDIF
      IF(F2*F1 .LE. F2*F3) THEN
        X3=X1
        F3=F1
      END IF
      X1=X
      F1=F2
 100  CONTINUE
C++ ENDE DER SCHLEIFE ++++++++++++++++++++++++++++++++++++++++++++++++
      IX=1
 999  RETURN
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION TPS5N(P,S)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision NULLS5N
      EXTERNAL NULLS5N
C
      EPS=1.D-6
C
      IF (P .GE. 5.D0) THEN
         IF (S .LT. 8.5D0) THEN
            T1=1073.15D0
            T2=1823.15D0
         ELSE
            T1=1323.15D0
            T2=2273.15D0
         ENDIF
      ELSE IF (P .GE. 1.D0) THEN
         IF (S .LT. 9.D0) THEN
            T1=1073.15D0
            T2=1823.15D0
         ELSE
            T1=1323.15D0
            T2=2273.15D0
         ENDIF
      ELSE IF (P .GE. 0.1D0) THEN
         IF (S .LT. 9.5D0) THEN
            T1=1073.15D0
            T2=1823.15D0
         ELSE
            T1=1073.15D0
            T2=2273.15D0
         ENDIF
      ELSE
         T1=1073.15D0
         T2=2273.15D0
      ENDIF
C
      CALL WNPS3(T1,T2,NULLS5N,S,P,EPS,X,IX)
C      
c      IF ((X .LE. 2273.15D0) .AND. (IX .EQ. 0)) THEN
      TPS5N=X
c      ELSE
c         TPS5N=0.D0
c      ENDIF
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NULLS5N(T,P,S)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      NULLS5N=SPT5N(P,T)-S
C
      END
C
C***********************************************************************
      SUBROUTINE REGSPV(P,V,IR3,IREG)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      COMMON/CGR/PGR,TGR13
C
      double precision POLD, VOLD
      SAVE POLD, VOLD
      INTEGER IR3OLD, IREGOLD
      SAVE IR3OLD, IREGOLD
      data POLD, VOLD, IR3OLD, IREGOLD / -1.D+0, -1.D+0, -1, -1 /
C
      IF ((DABS(V-VOLD).LT. 1.D-12) .AND.
     *   (DABS(P-POLD).LT. 1.D-12) .AND. (IR3 .EQ. IR3OLD)) THEN
           IREG=IREGOLD
         GOTO 999
      END IF
C          
      IREG=0
C
      IF ((P .GT. 100.D0) .OR. (P .LT. 5.D-4)) THEN
           IREG=0
        ELSEIF (P .GT. 22.064D0) THEN
           V1=VPT1N(P,273.15D0)
           V2=VPT1N(P,623.15D0)
           TG=FB23P(P)
           V3=VPT2N(P,TG)
           V4=VPT2N(P,1073.15D0)
           IF (V .LT. V1) THEN
              IREG=0
           ELSEIF (V .LT. V2) THEN
              IREG=1
           ELSEIF (V .LT. V3) THEN
              IREG=3
           ELSEIF (V .LT. V4) THEN
              IREG=2
           ELSE
              IREG=0
           ENDIF
        ELSEIF (P .GT. PGR) THEN
         TGR=TPVGR(P)
           V1=VPT1N(P,273.15D0)
           V2=VPT1N(P,623.15D0)
           TS=TSATPN(P)
           CALL FSATP(DV,DL,TOUT,TS,P)
           V3=1.D0/DL
           V4=1.D0/DV
           TG=FB23P(P)
           V5=VPT2N(P,TG)
           V6=VPT2N(P,1073.15D0)
         IF (TGR .LT. 273.15D0) THEN
              IF (V .LT. V1) THEN
                 IREG=0
              ELSEIF (V .LT. V2) THEN
                 IREG=1
              ELSEIF (V .LT. V3) THEN
                 IREG=3
              ELSEIF (V .LT. V4) THEN
                 IREG=9
              ELSEIF (V .LT. V5) THEN
                 IREG=3
              ELSEIF (V .LT. V6) THEN
                 IREG=2
              ELSE
                 IREG=0
              ENDIF
           ELSE
              VGR=VPT1N(P,TGR)
              IF (V .LT. VGR) THEN
                 IREG=0
              ELSEIF (V .LT. V2) THEN
                 IREG=1
              ELSEIF (V .LT. V3) THEN
                 IREG=3
              ELSEIF (V .LT. V4) THEN
                 IREG=9
              ELSEIF (V .LT. V5) THEN
                 IREG=3
              ELSEIF (V .LT. V6) THEN
                 IREG=2
              ELSE
                 IREG=0
              ENDIF
           ENDIF
        ELSEIF (P .GT. 10.D0) THEN
           TS=TSATPN(P)
           TGR=TPVGR(P)
           V2=VPT1N(P,TS)
           V3=VPT2N(P,TS)
           V4=VPT2N(P,1073.15D0)
           VGR=VPT1N(P,TGR)
           IF (V .LT. VGR) THEN
              IREG=0
           ELSEIF (V .LT. V2) THEN
              IREG=1
           ELSEIF (V .LT. V3) THEN
              IREG=9
           ELSEIF (V .LT. V4) THEN
              IREG=2
           ELSE
              IREG=0
           ENDIF
      ELSE
           TS=TSATPN(P)
           TGR=TPVGR(P)
         V1=VPT1N(P,273.15D0)
           V2=VPT1N(P,TS)
           V3=VPT2N(P,TS)
           V4=VPT2N(P,1073.15D0)
           V5=VPT5N(P,2273.15D0)
         IF (TGR .LE. TS) THEN
              VGR=VPT1N(P,TGR)
         ELSE
              VGR=V2
           ENDIF
           IF (V .LT. VGR) THEN
              IREG=0
           ELSEIF (V .LT. V2) THEN
              IREG=1
           ELSEIF (V .LT. V1) THEN
              IREG=1
           ELSEIF (V .LT. V3) THEN
              IREG=9
           ELSEIF (V .LT. V4) THEN
              IREG=2
           ELSEIF (V .LT. V5) THEN
              IREG=5
           ELSE
              IREG=0
           ENDIF
      ENDIF   
C
999   CONTINUE   
C
      POLD=P
      VOLD=V
          IR3OLD=IR3
      IREGOLD=IREG
C      
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION TPVGR(P)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
C      Gültig von 273.15 K bis 277.13 K
C
        a1 = -5.037450669   D-4
        a2 = -2.007841919       D-1
        a3 =  2.771333942       D2
C
        TPVGR = (a1*P+a2)*P+a3
C
        END
C
C***********************************************************************
      SUBROUTINE TPV(P,V,IREG,IP3,T1,T2)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
C oberer Wert
C
        IF (IP3 .EQ. 1) THEN
                IF (IREG .EQ. 0) THEN
                        CALL REGSPV(P,V,1,IREG)
                ENDIF
C  
                IF (IREG .EQ. 1) THEN
                        T1=TPV1N(P,V,1)
                ELSE IF (IREG .EQ. 2) THEN
                        T1=TPV2N(P,V)
                ELSE IF (IREG .EQ. 3) THEN
                        T1=TPV3N(P,V)
                ELSE IF (IREG .EQ. 5) THEN
                        T1=TPV5N(P,V)
                ELSE IF (IREG .EQ. 9) THEN
                        T1=TSATPN(P)
                ELSE
                        T1=0.D0
                ENDIF
                T2=0.D0
C
C unterer Wert
C
        ELSEIF (IP3 .EQ. 2) THEN
                IF (IREG .EQ. 0) THEN
                        CALL REGSPV(P,V,2,IREG)
                ENDIF
C  
                IF (IREG .EQ. 1) THEN
                        T1=TPV1N(P,V,2)
                ELSE IF (IREG .EQ. 2) THEN
                        T1=TPV2N(P,V)
                ELSE IF (IREG .EQ. 3) THEN
                        T1=TPV3N(P,V)
                ELSE IF (IREG .EQ. 5) THEN
                        T1=TPV5N(P,V)
                ELSE IF (IREG .EQ. 9) THEN
                        T1=TSATPN(P)
                ELSE
                        T1=0.D0
                ENDIF
                T2=0.D0
C
C beide Werte
C
        ELSE
                IF (IREG .EQ. 0) THEN
                        IRZ=0
                        CALL REGSPV(P,V,1,IREGU)
                ENDIF
C  
                IF (IREGU .EQ. 1) THEN
                        T1=TPV1N(P,V,1)
                ELSE IF (IREGU .EQ. 2) THEN
                        T1=TPV2N(P,V)
                ELSE IF (IREGU .EQ. 3) THEN
                        T1=TPV3N(P,V)
                ELSE IF (IREGU .EQ. 5) THEN
                        T1=TPV5N(P,V)
                ELSE IF (IREGU .EQ. 9) THEN
                        T1=TSATPN(P)
                ELSE
                        T1=0.D0
                ENDIF
C
                IF (IRZ .EQ. 0) THEN
                        CALL REGSPV(P,V,2,IREGO)
                ENDIF
C  
                IF (IREGO .EQ. 1) THEN
                        T2=TPV1N(P,V,2)
                ELSE
                        T2=0.D0
                ENDIF
C
        ENDIF
C
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION TPV1N(P,V,I3PE)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision NUPV1N
      EXTERNAL NUPV1N
C
      EPS=1.D-6
        TGR=TPVGR(P)
C
C oberer Wert
C
        IF (I3PE .EQ. 1) THEN
                IF ((TGR .GE. 273.15D0) .AND. (TGR .LE. 277.13D0)) THEN
                        T1=TGR
                        T2=TSATPN(P)
                ELSE
                        PGR=PSATTN(623.15D0)
                                IF (P .GT. PGR) THEN
                                        T1=273.15D0
                                        T2=623.15D0
                                ELSE
                                        T1=273.15D0
                                        T2=TSATPN(P)
                                ENDIF
                ENDIF
C
C unterer Wert
C
        ELSE
                IF ((TGR .GE. 273.15D0) .AND. (TGR .LE. 277.13D0)) THEN
                        T1=273.15D0
                        T2=TGR
                ELSE
                        PGR=PSATTN(623.15D0)
                                IF (P .GT. PGR) THEN
                                        T1=273.15D0
                                        T2=623.15D0
                                ELSE
                                        T1=273.15D0
                                        T2=TSATPN(P)
                                ENDIF
                ENDIF
        ENDIF
C
      CALL WNPV3(T1,T2,NUPV1N,V,P,EPS,X,IX)
C
      TPV1N=X
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUPV1N(T,P,V)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      NUPV1N=VPT1N(P,T)-V
C
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION TPV2N(P,V)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision NUPV2N
      EXTERNAL NUPV2N
      COMMON/CGR/PGR,TGR13
C
      EPS=1.D-6
C
      IF (P .GT. PGR) THEN
         T1=FB23P(P)
         T2=1073.15D0
      ELSE
         T1=TSATPN(P)
         T2=1073.15D0
      ENDIF
C
      CALL WNPV3(T1,T2,NUPV2N,V,P,EPS,X,IX)
C
      TPV2N=X
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUPV2N(T,P,V)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      NUPV2N=VPT2N(P,T)-V
C
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION TPV3N(P,V)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision NUPV3N
      EXTERNAL NUPV3N
C
      EPS=1.D-6
C
      T1=623.15D0
      T2=FB23P(P)
C
      CALL WNPV3(T1,T2,NUPV3N,V,P,EPS,X,IX)
C
      TPV3N=X
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUPV3N(T,P,V)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      NUPV3N=PVT3N(V,T)-P
C
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION TPV5N(P,V)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision NUPV5N
      EXTERNAL NUPV5N
C
      EPS=1.D-6
C
      T1=1073.15D0
      T2=2273.15D0
C
      CALL WNPV3(T1,T2,NUPV5N,V,P,EPS,X,IX)
C
      TPV5N=X
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUPV5N(T,P,V)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      NUPV5N=VPT5N(P,T)-V
C
      END
C
C*******************************************************************
      SUBROUTINE WNPV3(XA,XB,F,P,T,EPS,X,IX)
C*******************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      EXTERNAL F
      X1=XA
      F1=F(X1,T,P)
      X3=XB
      F3=F(X3,T,P)
      IX=0
C++ SCHLEIFE ZUR ERMITTLUNG DER NULLSTELLE +++++++++++++++++++++++++++
      DO 100 I=1,40
      IF( F1 .NE. F3 ) THEN
       X=X1+(X3-X1)*F1/(F1-F3)
      ELSE
       IX=3
       GOTO 999
      END IF
      IF(X .LT. 0.D0) X=(X1+X3)/2.D0
      IF(DABS(X) .LT. 1.D-8) THEN
        IF(DABS(X-X1).LT.EPS) RETURN
      ELSE
        IF(DABS((X-X1)/X).LT.EPS) RETURN
      END IF
      F2=F(X,T,P)
      X2=X1-(X1-X3)/2.D0
      P1=F2*F1
      P3=F2*F3
      IF ((P1 .LT. 0.D0) .AND. (P3 .LT. 0.D0)) THEN
        P1=DABS(P1)
        P3=DABS(P3)
      ENDIF
      IF(P1 .LE. P3) THEN
        X3=X1
        F3=F1
      END IF
      X1=X
      F1=F2
      IF((X2-X3)*(X2-X1) .GE. 0.D0) GOTO 100
c      X=(X1+X3)/2.D0
      IF ((DABS(F1/F3)) .LE. (0.3D0)) THEN
         X=X1+DABS(F1/F3)*1.5D0*(X3-X1)
         IREM=1
      ELSE IF ((DABS(F1/F3)) .GE. (3.D0)) THEN
         IREM=2
      ELSE
         X=(X1+X3)/2.D0
         IREM=0
      ENDIF
      F2=F(X,T,P)
      IF ((((F2*F1) .GE. 0.D0) .AND. (IREM .EQ. 1)) .OR.
     1    (((F2*F3) .GE. 0.D0) .AND. (IREM .EQ. 2))) THEN
         X=(X1+X3)/2.D0
         F2=F(X,T,P)
      ENDIF
      IF(F2*F1 .LE. F2*F3) THEN
        X3=X1
        F3=F1
      END IF
      X1=X
      F1=F2
 100  CONTINUE
C++ ENDE DER SCHLEIFE ++++++++++++++++++++++++++++++++++++++++++++++++
      IX=1
 999  RETURN
      END
C
C***********************************************************************
      SUBROUTINE REGSTH(T,H,IR1,IREG)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      double precision TOLD, HOLD
      INTEGER IR1OLD, IREGOLD
      SAVE TOLD, HOLD
      SAVE IR1OLD, IREGOLD
      data  TOLD, HOLD, IR1OLD, IREGOLD / -1.D+0, -1.D-6, -1, -1 /
C
      IF ((DABS(T-TOLD).LT. 1.D-6) .AND.
     *   (DABS(H-HOLD).LT. 1.D-6) .AND. (IR1 .EQ. IR1OLD)) THEN
         IREG=IREGOLD
         GOTO 999
      END IF
C          
      IF ((T .LT. 273.15D0) .OR. (T .GT. 2273.15D0)) THEN
         IREG=0
         GOTO 999
      ENDIF
C
        IREG=0
C 
      IF (T .LE. 623.15D0) THEN
         PG=PSATTN(T)
         H1=HPT1N(PG,T)
         H2=HPT2N(PG,T)
         IF ((H .LT. H2) .AND. (IR1 .EQ. 1)) THEN
            IF (T .LE. 520.767185D0) THEN
               H100=HPT1N(100.D0,T)
               IF ((H .GE. H1) .AND. (H .LE. H2)) THEN
                  IREG=9
               ELSEIF ((H .GE. H1) .AND. (H .LE. H100)) THEN
                  IREG=1
               ELSE
                  IREG=0
               ENDIF
ccc
            ELSEIF (T .LE. 611.9697430625D0) THEN
               IF ((H .GE. H1) .AND. (H .LE. H2)) THEN
                  IREG=9
               ELSE
                            PGH=PTHGR(T)
                    HMIN=HPT1N(PGH,T)
                            IF (H .GE. HMIN) THEN
                     IREG=1
                  ELSE
                     IREG=0
                  ENDIF
                 ENDIF
            ELSE
               IF ((H .GE. H1) .AND. (H .LE. H2)) THEN
                  IREG=9
               ELSE
                    HMIN=HPT1N(100.D0,T)
                  IF ((H .GE. HMIN) .AND. (H .LE. H1)) THEN
                     IREG=1
                  ELSE
                     IREG=0
                  ENDIF
               ENDIF
            ENDIF
         ELSEIF ((H .LT. H2) .AND. (IR1 .EQ. 2)) THEN
            IF ((H .GE. H1) .AND. (H .LE. H2) .AND.
     &       (H .GT. 2000.D0)) THEN
               IREG=9
            ELSEIF (T .LE. 520.767185D0) THEN
               H100=HPT1N(100.D0,T)
               IF ((H .GE. H1) .AND. (H .LE. H100) .AND. 
     &          (H .LE. 2000.D0)) THEN
                  PN=PTH1N(T,H,2)
                  IF (PN .LE. 102.D0) THEN
                     IREG=1
                  ELSE
                     IREG=0
                  ENDIF
               ELSE
                  IREG=0
               ENDIF
            ELSEIF (T .LE. 611.9697430625D0) THEN
                         PGH=PTHGR(T)
                 HMIN=HPT1N(PGH,T)
                 HMAX=HPT1N(100.D0,T)
                         IF ((H .GE. HMIN) .AND. (H .LE. HMAX)) THEN
                  IREG=1
               ELSE
                  IREG=0
               ENDIF
            ELSE
               IREG=0
            ENDIF
         ELSEIF (H .GE. H2) THEN
            HG=HPT2N(5.D-4,T)
            IF (H .LE. HG) THEN
               IREG=2
            ELSE
               IREG=0
            ENDIF
         ELSE
            IREG=9
         ENDIF
      ELSEIF (T .LE. 863.15D0) THEN
         PG=FB23(T)
         HB=HPT2N(PG,T)
         IF (H .GE. HB) THEN
            HG=HPT2N(5.D-4,T)
            IF (H .LE. HG) THEN
               IREG=2
            ELSE
               IREG=0
            ENDIF
         ELSE
            V100=VPT3N(100.D0,T)
            HG=HVT3N(V100,T)
            IF (H .LT. HG) THEN
               IREG=0
            ELSEIF (T .GT. 647.096D0) THEN
               IREG=3
            ELSE
               PG=PSATTN(T)
               CALL FSATP(DV,DL,TOUT,T,PG)
               V1=1.D0/DL
               V2=1.D0/DV
               H1=HVT3N(V1,T)
               H2=HVT3N(V2,T)
               IF ((H .GT. H1) .AND. (H .LT. H2)) THEN
                  IREG=9
               ELSE
                  IREG=3
               ENDIF
            ENDIF
         ENDIF
      ELSEIF (T .LE. 1073.15D0) THEN
         HO=HPT2N(100.D0,T)
         HU=HPT2N(5.D-4,T)
         IF ((H .GE. HO) .AND. (H .LE. HU)) THEN
            IREG=2
         ELSE
            IREG=0
         ENDIF
      ELSE
         HO=HPT5N(10.D0,T)
         HU=HPT5N(5.D-4,T)
         IF ((H .GE. HO) .AND. (H .LE. HU)) THEN
            IREG=5
         ELSE
            IREG=0
         ENDIF
      ENDIF
C
999   CONTINUE
C
      TOLD=T
      HOLD=H
        IR1OLD=IR1
      IREGOLD=IREG
C   
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION PTHGR(T)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
c      Gültig von 520.767185D0 K bis 611.9697430625D0 K
C
        a1=4.22963120715754d-08
        a2=-9.55898521095589d-05
        a3=0.080230855541771d0
        a4=-28.5874899158642d0
        a5=3522.33800777614d0
C
        pthgr= (((a1*t+a2)*t+a3)*t+a4)*t+a5
C
        END
C***********************************************************************
      SUBROUTINE PTH(T,H,IREG,IP1,P1,P2)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      COMMON/CCRIT/VC,SC,HC
C
      IF (IP1 .EQ. 1) THEN
C
ccc unterer Druck
C     
         IF (IREG .EQ. 0) THEN
            CALL REGSTH(T,H,1,IREG)
         ENDIF
C  
         IF (IREG .EQ. 1) THEN
            P1=PTH1N(T,H,1)
         ELSE IF (IREG .EQ. 2) THEN
            P1=PTH2N(T,H)
         ELSE IF (IREG .EQ. 3) THEN
            P1=PTH3N(T,H)
         ELSE IF (IREG .EQ. 5) THEN
            P1=PTH5N(T,H)
         ELSE IF (IREG .EQ. 9) THEN
            P1=PSATTN(T)
           ELSE
                  P1=0.D0
         ENDIF
         P2=0.D0
C
      ELSEIF (IP1 .EQ. 2) THEN
C
ccc oberer Druck
C     
         IF (IREG .EQ. 0) THEN
            CALL REGSTH(T,H,2,IREG)
         ENDIF
C  
         IF (IREG .EQ. 1) THEN
            P1=PTH1N(T,H,2)
         ELSE IF (IREG .EQ. 2) THEN
            P1=PTH2N(T,H)
         ELSE IF (IREG .EQ. 3) THEN
            P1=PTH3N(T,H)
         ELSE IF (IREG .EQ. 5) THEN
            P1=PTH5N(T,H)
         ELSE IF ((IREG .EQ. 9) .AND. (H .GE. HC)) THEN
            P1=PSATTN(T)
           ELSE
                  P1=0.D0
         ENDIF
         P2=0.D0
C
      ELSE
C
ccc beide Drücke
C     
         IF (IREG .EQ. 0) THEN
            IRZ=0
            CALL REGSTH(T,H,1,IREGU)
         ENDIF
C  
         IF (IREGU .EQ. 1) THEN
            P1=PTH1N(T,H,1)
         ELSE IF (IREGU .EQ. 2) THEN
            P1=PTH2N(T,H)
         ELSE IF (IREGU .EQ. 3) THEN
            P1=PTH3N(T,H)
         ELSE IF (IREGU .EQ. 5) THEN
            P1=PTH5N(T,H)
         ELSE IF (IREGU .EQ. 9) THEN
            P1=PSATTN(T)
         ELSE
            P1=0.D0
         ENDIF   
         IF (IRZ .EQ. 0) THEN   
            CALL REGSTH(T,H,2,IREGO)
         ENDIF
         IF (IREGO .EQ. 1) THEN
            P2=PTH1N(T,H,2)
         ELSE IF (IREGO .EQ. 2) THEN
            P2=0.D0
         ELSE IF (IREGO .EQ. 3) THEN
            P2=0.d0
         ELSE IF (IREGO .EQ. 5) THEN
            P2=0.d0
         ELSE IF (IREGO .EQ. 9) THEN
            P2=0.d0
         ELSE 
            P2=0.d0
         ENDIF   
C         
      ENDIF   
C
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION PTH1N(T,H,I1PE)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision NUTH1U,NUTH1O
      EXTERNAL NUTH1U,NUTH1O
C
      EPS=1.D-9
C
      IF (I1PE .EQ. 1) THEN
         IF (T .LE. 520.767185D0) THEN
            HZ1=HPT1N(60.D0,T)
            HZ2=HPT1N(25.D0,T)
            IF (H .GE. HZ1) THEN
               P1=60.D0
               P2=100.D0
            ELSE IF (H .GE. HZ2) THEN
               P1=25.D0
               P2=60.D0
            ELSE
               P1=PSATTN(T)
               P2=25.D0
            ENDIF
         ELSEIF (T .LE. 611.9697430625D0) THEN
              P1=PSATTN(T)
              P2=PTHGR(T)
         ELSE
            HZ1=HPT1N(50.D0,T)
            IF (H .LE. HZ1) THEN
               P1=50.D0
               P2=80.D0
            ELSE
               P1=PSATTN(T)
               P2=50.D0
            ENDIF
         ENDIF
ccc      
         CALL WNTH3(P1,P2,NUTH1U,H,T,EPS,X,IX)
C
      ELSE
         IF (T .LE. 520.767185D0) THEN
            HZ1=HPT1N(60.D0,T)
            HZ2=HPT1N(25.D0,T)
            IF (H .GE. HZ1) THEN
               P1=60.D0
               P2=100.D0
            ELSE IF (H .GE. HZ2) THEN
               P1=25.D0
               P2=60.D0
            ELSE
               P1=PSATTN(T)
               P2=25.D0
            ENDIF
         ELSEIF (T .LE. 611.9697430625D0) THEN
              P1=PTHGR(T)
              P2=100.D0
         ELSE
            HZ1=HPT1N(60.D0,T)
            IF (H .LE. HZ1) THEN
               P1=95.D0
               P2=100.D0
            ELSE
               P1=PSATTN(T)
               P2=60.D0
            ENDIF
         ENDIF
C
         CALL WNTH3(P1,P2,NUTH1O,H,T,EPS,X,IX)
C
      ENDIF
C
      PTH1N=X
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUTH1U(P,T,H)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      NUTH1U=HPT1N(P,T)-H
C
      END
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUTH1O(P,T,H)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      NUTH1O=HPT1N(P,T)-H
C
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION PTH2N(T,H)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision NUTH2N
      EXTERNAL NUTH2N
C
      EPS=1.D-6
C
      IF (T .LE. 623.15D0) THEN
         P1=5.D-4
         P2=PSATTN(T)
      ELSEIF (T .LE. 863.15D0) THEN
         HZ=HPT2N(15.D0,T)
         IF (H .GE. HZ) THEN
            P1=5.D-4
            P2=15.D0
         ELSE
            P1=15.D0
            P2=FB23(T)
         ENDIF
      ELSE
         HZ1=HPT2N(60.D0,T)
         HZ2=HPT2N(20.D0,T)
         IF (H .LE. HZ1) THEN
            P1=60.D0
            P2=100.D0
         ELSEIF (H .LE. HZ2) THEN
            P1=20.D0
            P2=60.D0
         ELSE
            P1=5.D-4
            P2=20.D0
         ENDIF
      ENDIF   
C
      CALL WNTH3(P1,P2,NUTH2N,H,T,EPS,X,IX)
C
      PTH2N=X
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUTH2N(P,T,H)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      NUTH2N=HPT2N(P,T)-H
C
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION PTH3N(T,H)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      COMMON/CSUB2/TC,PC,DC
      double precision NUTH3N
      EXTERNAL NUTH3N
C
      EPS=1.D-6
C
      IF (T .LE. TC) THEN
         V60=VPT3N(60.D0,T)
         HZ1=HVT3N(V60,T)
         PG=PSATTN(T)
         CALL FSATP(DV,DL,TOUT,T,PG)
               V1=1.D0/DL
         HZ2=HVT3N(V1,T)
         IF (H .LE. HZ1) THEN
            P1=60.D0
            P2=100.D0
         ELSE IF (H .LE. HZ2) THEN
            P1=PG+0.01D0
            P2=60.D0
         ELSE
            P1=FB23(T)
            P2=PG-0.01D0
         ENDIF
      ELSE
         P1=FB23(T)
         P2=100.D0
      ENDIF
C
      CALL WNTH3(P1,P2,NUTH3N,H,T,EPS,X,IX)
C
      PTH3N=X
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUTH3N(P,T,H)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      V=VPT3N(P,T)
      NUTH3N=HVT3N(V,T)-H
C
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION PTH5N(T,H)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision NUTH5N
      EXTERNAL NUTH5N
C
      EPS=1.D-6
C
      HZ1=HPT5N(3.D0,T)
      HZ2=HPT5N(0.1D0,T)
      IF (H .LE. HZ1) THEN
         P1=3.D0
         P2=10.D0
      ELSE IF (H .LE. HZ2) THEN
         P1=0.1D0
         P2=3.D0
      ELSE
         P1=5.D-4
         P2=0.1D0
      ENDIF
C
      CALL WNTH3(P1,P2,NUTH5N,H,T,EPS,X,IX)
C
      PTH5N=X
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUTH5N(P,T,H)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      NUTH5N=HPT5N(P,T)-H
C
      END
C
C*******************************************************************
      SUBROUTINE WNTH3(XA,XB,F,P,T,EPS,X,IX)
C*******************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      EXTERNAL F
      X1=XA
      F1=F(X1,T,P)
      X3=XB
      F3=F(X3,T,P)
      IX=0
C++ SCHLEIFE ZUR ERMITTLUNG DER NULLSTELLE +++++++++++++++++++++++++++
      DO 100 I=1,40
      IF( F1 .NE. F3 ) THEN
       X=X1+(X3-X1)*F1/(F1-F3)
      ELSE
       IX=3
       GOTO 999
      END IF
      IF(X .LT. 0.D0) X=(X1+X3)/2.D0
      IF(DABS(X) .LT. 1.D-8) THEN
        IF(DABS(X-X1).LT.EPS) RETURN
      ELSE
        IF(DABS((X-X1)/X).LT.EPS) RETURN
      END IF
      F2=F(X,T,P)
      X2=X1-(X1-X3)/2.D0
      P1=F2*F1
      P3=F2*F3
      IF ((P1 .LT. 0.D0) .AND. (P3 .LT. 0.D0)) THEN
        P1=DABS(P1)
        P3=DABS(P3)
      ENDIF
      IF(P1 .LE. P3) THEN
        X3=X1
        F3=F1
      END IF
      X1=X
      F1=F2
      IF((X2-X3)*(X2-X1) .GE. 0.D0) GOTO 100
c      X=(X1+X3)/2.D0
      IF ((DABS(F1/F3)) .LE. (0.3D0)) THEN
         X=X1+DABS(F1/F3)*1.5D0*(X3-X1)
         IREM=1
      ELSE IF ((DABS(F1/F3)) .GE. (3.D0)) THEN
         IREM=2
      ELSE
         X=(X1+X3)/2.D0
         IREM=0
      ENDIF
      F2=F(X,T,P)
      IF ((((F2*F1) .GE. 0.D0) .AND. (IREM .EQ. 1)) .OR.
     1    (((F2*F3) .GE. 0.D0) .AND. (IREM .EQ. 2))) THEN
         X=(X1+X3)/2.D0
         F2=F(X,T,P)
      ENDIF
      IF(F2*F1 .LE. F2*F3) THEN
        X3=X1
        F3=F1
      END IF
      X1=X
      F1=F2
 100  CONTINUE
C++ ENDE DER SCHLEIFE ++++++++++++++++++++++++++++++++++++++++++++++++
      IX=1
 999  RETURN
      END
C
C***********************************************************************
      SUBROUTINE REGSTS(T,S,IR2,IREG)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      double precision TOLD, SOLD
      SAVE TOLD, SOLD
      INTEGER IR2OLD, IREGOLD
      SAVE IR2OLD, IREGOLD
      data TOLD, SOLD, IR2OLD, IREGOLD / -1.D+0, -1.D-6, -1, -1 /
C
      IF ((DABS(T-TOLD).LT. 1.D-6) .AND.
     *   (DABS(S-SOLD).LT. 1.D-6) .AND. (IR2 .EQ. IR2OLD)) THEN
         IREG=IREGOLD
         GOTO 999
      END IF
C          
      IF ((T .LT. 273.15D0) .OR. (T .GT. 2273.15D0)) THEN
         IREG=0
         GOTO 999
      ENDIF
C
      IREG=0
C
      IF (T .LE. 623.15D0) THEN
        PG=PSATTN(T)
        S1=SPT1N(PG,T)
        S2=SPT2N(PG,T)
        IF (T .LE. 277.13D0) THEN
C       
C oberer Wert
C
            IF ((IR2 .EQ. 1) .AND. (S .LE. S2)) THEN
                S100=SPT1N(100.D0,T)
                PSG=PTSGR(T)
                SG=SPT1N(PSG,T)
                IF ((S .GE. S100) .AND. (S .LE. SG)) THEN
                    IREG=1
                ELSEIF ((S .GE. S1) .AND. (S .LE. S2)) THEN
                    IREG=9
                ELSE
                    IREG=0
                ENDIF                   
C
C unterer Wert
C
            ELSEIF ((IR2 .EQ. 2) .AND. (S .LE. S2)) THEN
                PSG=PTSGR(T)
                SG=SPT1N(PSG,T)
                IF ((S .LE. SG) .AND. (S .GT. S1)) THEN
                    IREG=1
                ELSEIF ((S .GE. S1) .AND. (S .LE. S2)) THEN
                    IREG=9
                ELSE
                    IREG=0
                ENDIF
            ELSEIF (S .GT. S2) THEN
                SG=SPT2N(5.D-4,T)
                IF (S .LE. SG) THEN
                    IREG=2
                ELSE
                    IREG=0
                ENDIF
            ELSE
                IREG=0                          
            ENDIF
        ELSE
            IF (S .LE. S1) THEN
                SG=SPT1N(100.D0,T)
                IF (S .GE. SG) THEN
                    IREG=1
                ELSE
                    IREG=0
                ENDIF
            ELSE
                S2=SPT2N(PG,T)
                IF (S .GE. S2) THEN
                    SG=SPT2N(5.D-4,T)
                    IF (S .LE. SG) THEN
                        IREG=2
                    ELSE
                        IREG=0
                    ENDIF
                ELSE
                        IREG=9
                ENDIF
            ENDIF
        ENDIF
C
      ELSEIF (T .LE. 863.15D0) THEN
         PG=FB23(T)
         SB=SPT2N(PG,T)
         IF (S .GE. SB) THEN
            SG=SPT2N(5.D-4,T)
            IF (S .LE. SG) THEN
               IREG=2
            ELSE
               IREG=0
            ENDIF
         ELSE
            V100=VPT3N(100.D0,T)
            SG=SVT3N(V100,T)
            IF (S .LT. SG) THEN
               IREG=0
            ELSEIF (T .GT. 647.096D0) THEN
               IREG=3
            ELSE
               PG=PSATTN(T)
               CALL FSATP(DV,DL,TOUT,T,PG)
               V1=1.D0/DL
               V2=1.D0/DV
               S1=SVT3N(V1,T)
               S2=SVT3N(V2,T)
               IF ((S .GT. S1) .AND. (S .LT. S2)) THEN
                  IREG=9
               ELSE
                  IREG=3
               ENDIF
            ENDIF
         ENDIF
      ELSEIF (T .LE. 1073.15D0) THEN
         SO=SPT2N(100.D0,T)
         SU=SPT2N(5.D-4,T)
         IF ((S .GE. SO) .AND. (S .LE. SU)) THEN
            IREG=2
         ELSE
            IREG=0
         ENDIF
      ELSE
         SO=SPT5N(10.D0,T)
         SU=SPT5N(5.D-4,T)
         IF ((S .GE. SO) .AND. (S .LE. SU)) THEN
            IREG=5
         ELSE
            IREG=0
         ENDIF
      ENDIF
C
999   CONTINUE   
C
      TOLD=T
      SOLD=S
          IR2OLD=IR2
      IREGOLD=IREG
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION PTSGR(T)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
c      Gültig von 273.15 K bis 277.13 K
C
        a1 = -2.429465065   D-05
        a2 =  0.02563798356 D0
        a3 = -1.018127106   D01
        a4 =  1.799358124   D03
        a5 = -1.190999198   D05
C
        PTSGR = (((a1*T+a2)*T+a3)*T+a4)*T+a5
C
        END
C
C***********************************************************************
      SUBROUTINE PTS(T,S,IREG,IP2,P1,P2)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
C      
C oberer Druck
C
        IF (IP2 .EQ. 1) THEN
                IF (IREG .EQ. 0) THEN
                        CALL REGSTS(T,S,1,IREG)
                ENDIF
C  
                IF (IREG .EQ. 1) THEN
                        P1=PTS1N(T,S,1)
                ELSE IF (IREG .EQ. 2) THEN
                        P1=PTS2N(T,S)
                ELSE IF (IREG .EQ. 3) THEN
                        P1=PTS3N(T,S)
                ELSE IF (IREG .EQ. 5) THEN
                        P1=PTS5N(T,S)
                ELSE IF (IREG .EQ. 9) THEN
                        P1=PSATTN(T)
                ELSE
                        P1=0.D0
                ENDIF
                P2=0.D0
C      
C unterer Druck
C
        ELSEIF (IP2 .EQ. 2) THEN
                IF (IREG .EQ. 0) THEN
                        CALL REGSTS(T,S,2,IREG)
                ENDIF
C  
                IF (IREG .EQ. 1) THEN
                        P1=PTS1N(T,S,2)
                ELSE IF (IREG .EQ. 2) THEN
                        P1=PTS2N(T,S)
                ELSE IF (IREG .EQ. 3) THEN
                        P1=PTS3N(T,S)
                ELSE IF (IREG .EQ. 5) THEN
                        P1=PTS5N(T,S)
                ELSE IF (IREG .EQ. 9) THEN
                        P1=PSATTN(T)
                ELSE
                        P1=0.D0
                ENDIF
                P2=0.D0
C      
C beide Drücke
C
        ELSE
                IF (IREG .EQ. 0) THEN
                        IRZ=0
                        CALL REGSTS(T,S,1,IREGU)
                ENDIF
C  
                IF (IREGU .EQ. 1) THEN
                        P1=PTS1N(T,S,1)
                ELSE IF (IREGU .EQ. 2) THEN
                        P1=PTS2N(T,S)
                ELSE IF (IREGU .EQ. 3) THEN
                        P1=PTS3N(T,S)
                ELSE IF (IREGU .EQ. 5) THEN
                        P1=PTS5N(T,S)
                ELSE IF (IREG .EQ. 9) THEN
                        P1=PSATTN(T)
                ELSE
                        P1=0.D0
                ENDIF
C
                IF (IRZ .EQ. 0) THEN
                        CALL REGSTS(T,S,2,IREGO)
                ENDIF
C  
                IF (IREGO .EQ. 1) THEN
                        P2=PTS1N(T,S,2)
                ELSE
                        P2=0.D0
                ENDIF
C      
        ENDIF
C        
        END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION PTS1N(T,S,I2PE)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision NUTS1N
      EXTERNAL NUTS1N
C
      EPS=1.D-6
C
C oberer Wert
C       
        IF (I2PE .EQ. 1) THEN
                IF (T .LE. 277.13D0) THEN
                        P1=PTSGR(T)
                        P2=100.D0
                ELSE
                        SZ1=SPT1N(60.D0,T)
                        SZ2=SPT1N(25.D0,T)
                        IF (S .LE. SZ1) THEN
                                P1=60.D0
                                P2=100.D0
                        ELSE IF (S .LE. SZ2) THEN
                                P1=25.D0
                                P2=60.D0
                        ELSE
                                P1=PSATTN(T)
                                P2=25.D0
                        ENDIF
                ENDIF
                CALL WNTS3(P1,P2,NUTS1N,S,T,EPS,X,IX)
C
C unterer Wert
C
        ELSE
                IF (T .LE. 277.13D0) THEN
                        P1=PSATTN(T)
                        P2=PTSGR(T)
                ELSE
                        SZ1=SPT1N(60.D0,T)
                        SZ2=SPT1N(25.D0,T)
                        IF (S .LE. SZ1) THEN
                                P1=60.D0
                                P2=100.D0
                        ELSE IF (S .LE. SZ2) THEN
                                P1=25.D0
                                P2=60.D0
                        ELSE
                                P1=PSATTN(T)
                                P2=25.D0
                        ENDIF
                ENDIF
                CALL WNTS3(P1,P2,NUTS1N,S,T,EPS,X,IX)
        ENDIF   
C
C
      PTS1N=X
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUTS1N(P,T,S)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      NUTS1N=SPT1N(P,T)-S
C
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION PTS2N(T,S)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision NUTS2N
      EXTERNAL NUTS2N
C
      EPS=1.D-6
C
      IF (T .LE. 623.15D0) THEN
         P1=5.D-4
         P2=PSATTN(T)
      ELSEIF (T .LE. 863.15D0) THEN
         SZ=SPT2N(15.D0,T)
         IF (S .GE. SZ) THEN
            P1=5.D-4
            P2=15.D0
         ELSE
            P1=15.D0
            P2=FB23(T)
         ENDIF
      ELSE
         SZ1=SPT2N(60.D0,T)
         SZ2=SPT2N(20.D0,T)
         IF (S .LE. SZ1) THEN
            P1=60.D0
            P2=100.D0
         ELSEIF (S .LE. SZ2) THEN
            P1=20.D0
            P2=60.D0
         ELSE
            P1=5.D-4
            P2=20.D0
         ENDIF
      ENDIF   
C
      CALL WNTS3(P1,P2,NUTS2N,S,T,EPS,X,IX)
C
      PTS2N=X
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUTS2N(P,T,S)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      NUTS2N=SPT2N(P,T)-S
C
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION PTS3N(T,S)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      COMMON/CSUB2/TC,PC,DC
      double precision NUTS3N
      EXTERNAL NUTS3N
C
      EPS=1.D-6
C
      IF (T .LE. TC) THEN
         V60=VPT3N(60.D0,T)
         SZ1=SVT3N(V60,T)
         PG=PSATTN(T)
         CALL FSATP(DV,DL,TOUT,T,PG)
               V1=1.D0/DL
         SZ2=SVT3N(V1,T)
         IF (S .LE. SZ1) THEN
            P1=60.D0
            P2=100.D0
         ELSE IF (S .LE. SZ2) THEN
            P1=PG+0.01D0
            P2=60.D0
         ELSE
            P1=FB23(T)
            P2=PG-0.01D0
         ENDIF
      ELSE
         P1=FB23(T)
         P2=100.D0
      ENDIF
C
      CALL WNTS3(P1,P2,NUTS3N,S,T,EPS,X,IX)
C
      PTS3N=X
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUTS3N(P,T,S)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      V=VPT3N(P,T)
      NUTS3N=SVT3N(V,T)-S
C
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION PTS5N(T,S)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision NUTS5N
      EXTERNAL NUTS5N
C
      EPS=1.D-6
C
      SZ1=SPT5N(3.D0,T)
      SZ2=SPT5N(0.1D0,T)
      IF (S .LE. SZ1) THEN
         P1=3.D0
         P2=10.D0
      ELSE IF (S .LE. SZ2) THEN
         P1=0.1D0
         P2=3.D0
      ELSE
         P1=5.D-4
         P2=0.1D0
      ENDIF
C
      CALL WNTS3(P1,P2,NUTS5N,S,T,EPS,X,IX)
C
      PTS5N=X
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUTS5N(P,T,S)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      NUTS5N=SPT5N(P,T)-S
C
      END
C
C*******************************************************************
      SUBROUTINE WNTS3(XA,XB,F,P,T,EPS,X,IX)
C*******************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      EXTERNAL F
      X1=XA
      F1=F(X1,T,P)
      X3=XB
      F3=F(X3,T,P)
      IX=0
C++ SCHLEIFE ZUR ERMITTLUNG DER NULLSTELLE +++++++++++++++++++++++++++
      DO 100 I=1,40
      IF( F1 .NE. F3 ) THEN
       X=X1+(X3-X1)*F1/(F1-F3)
      ELSE
       IX=3
       GOTO 999
      END IF
      IF(X .LT. 0.D0) X=(X1+X3)/2.D0
      IF(DABS(X) .LT. 1.D-8) THEN
        IF(DABS(X-X1).LT.EPS) RETURN
      ELSE
        IF(DABS((X-X1)/X).LT.EPS) RETURN
      END IF
      F2=F(X,T,P)
      X2=X1-(X1-X3)/2.D0
      P1=F2*F1
      P3=F2*F3
      IF ((P1 .LT. 0.D0) .AND. (P3 .LT. 0.D0)) THEN
        P1=DABS(P1)
        P3=DABS(P3)
      ENDIF
      IF(P1 .LE. P3) THEN
        X3=X1
        F3=F1
      END IF
      X1=X
      F1=F2
      IF((X2-X3)*(X2-X1) .GE. 0.D0) GOTO 100
c      X=(X1+X3)/2.D0
      IF ((DABS(F1/F3)) .LE. (0.3D0)) THEN
         X=X1+DABS(F1/F3)*1.5D0*(X3-X1)
         IREM=1
      ELSE IF ((DABS(F1/F3)) .GE. (3.D0)) THEN
         IREM=2
      ELSE
         X=(X1+X3)/2.D0
         IREM=0
      ENDIF
      F2=F(X,T,P)
      IF ((((F2*F1) .GE. 0.D0) .AND. (IREM .EQ. 1)) .OR.
     1    (((F2*F3) .GE. 0.D0) .AND. (IREM .EQ. 2))) THEN
         X=(X1+X3)/2.D0
         F2=F(X,T,P)
      ENDIF
      IF(F2*F1 .LE. F2*F3) THEN
        X3=X1
        F3=F1
      END IF
      X1=X
      F1=F2
 100  CONTINUE
C++ ENDE DER SCHLEIFE ++++++++++++++++++++++++++++++++++++++++++++++++
      IX=1
 999  RETURN
      END
C
C***********************************************************************
      SUBROUTINE REGSTV(T,V,IREG)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      double precision TOLD, VOLD
      SAVE TOLD, VOLD
      INTEGER IREGOLD
      SAVE IREGOLD
      data TOLD, VOLD, IREGOLD / -1.D+0, -1.D+0, -1 /
C
      IF ((DABS(T-TOLD).LT. 1.D-6) .AND.
     *   (DABS(V-VOLD).LT. 1.D-6)) THEN
         IREG=IREGOLD
         GOTO 999
      END IF
C          
      IF ((T .LT. 273.15D0) .OR. (T .GT. 2273.15D0)) THEN
         IREG=0
         GOTO 999
      ENDIF
C
        IREG=0
C
      IF (T .LE. 623.15D0) THEN
         PG=PSATTN(T)
         V1=VPT1N(PG,T)
         IF (V .LE. V1) THEN
            VG=VPT1N(101.D0,T)
            IF (V .GE. VG) THEN
               IREG=1
            ELSE
               IREG=0
            ENDIF
         ELSE
            V2=VPT2N(PG,T)
            IF (V .GE. V2) THEN
               VG=VPT2N(5.D-4,T)
               IF (V .LE. VG) THEN
                  IREG=2
               ELSE
                  IREG=0
               ENDIF
            ELSE
               IREG=9
            ENDIF
         ENDIF
      ELSEIF (T .LE. 863.15D0) THEN
         PG=FB23(T)
         VB=VPT2N(PG,T)
         IF (V .GE. VB) THEN
            VG=VPT2N(5.D-4,T)
            IF (V .LE. VG) THEN
               IREG=2
            ELSE
               IREG=0
            ENDIF
         ELSE
            VG=VPT3N(101.D0,T)
            IF (V .LE. VG) THEN
               IREG=0
            ELSEIF (T .GT. 647.096D0) THEN
               IREG=3
            ELSE
               PG=PSATTN(T)
               CALL FSATP(DV,DL,TOUT,T,PG)
               V1=1.D0/DL
               V2=1.D0/DV
               IF ((V .GT. V1) .AND. (V .LT. V2)) THEN
                  IREG=9
               ELSE
                  IREG=3
               ENDIF
            ENDIF
         ENDIF
      ELSEIF (T .LE. 1073.15D0) THEN
         VO=VPT2N(101.D0,T)
         VU=VPT2N(5.D-4,T)
         IF ((V .GE. VO) .AND. (V .LE. VU)) THEN
            IREG=2
         ELSE
            IREG=0
         ENDIF
      ELSE
         VO=VPT5N(10.1D0,T)
         VU=VPT5N(5.D-4,T)
         IF ((V .GE. VO) .AND. (V .LE. VU)) THEN
            IREG=5
         ELSE
            IREG=0
         ENDIF
      ENDIF
C
999   CONTINUE 
C
      TOLD=T
      VOLD=V
      IREGOLD=IREG
C         
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION PTV(T,V,IREG)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      IF (IREG .EQ. 0) THEN
         CALL REGSTV(T,V,IREG)
      ENDIF
C  
      IF (IREG .EQ. 1) THEN
         PTV=PTV1N(T,V)
      ELSE IF (IREG .EQ. 2) THEN
         PTV=PTV2N(T,V)
      ELSE IF (IREG .EQ. 3) THEN
         PTV=PVT3N(V,T)
      ELSE IF (IREG .EQ. 5) THEN
         PTV=PTV5N(T,V)
      ELSE IF (IREG .EQ. 9) THEN
         PTV=PSATTN(T)
      ELSE
         PTV=0.D0
      ENDIF
C
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION PTV1N(T,V)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision NUTV1N
      EXTERNAL NUTV1N
C
      EPS=1.D-6
C
      VZ1=VPT1N(60.D0,T)
      VZ2=VPT1N(25.D0,T)
      IF (V .LE. VZ1) THEN
         P1=60.D0
         P2=100.D0
      ELSE IF (V .LE. VZ2) THEN
         P1=25.D0
         P2=60.D0
      ELSE
         P1=PSATTN(T)
         P2=25.D0
      ENDIF
C
      CALL WNTV3(P1,P2,NUTV1N,V,T,EPS,X,IX)
C
      PTV1N=X
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUTV1N(P,T,V)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      NUTV1N=VPT1N(P,T)-V
C
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION PTV2N(T,V)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision NUTV2N
      EXTERNAL NUTV2N
C
      EPS=1.D-9
C
      IF (T .LE. 623.15D0) THEN
         P1=5.D-4
         P2=PSATTN(T)
      ELSEIF (T .LE. 863.15D0) THEN
         VZ=VPT2N(15.D0,T)
         IF (V .GE. VZ) THEN
            P1=5.D-4
            P2=15.D0
         ELSE
            P1=15.D0
            P2=FB23(T)
         ENDIF
      ELSE
         VZ1=VPT2N(60.D0,T)
         VZ2=VPT2N(20.D0,T)
         IF (V .LE. VZ1) THEN
            P1=60.D0
            P2=100.D0
         ELSEIF (V .LE. VZ2) THEN
            P1=20.D0
            P2=60.D0
         ELSE
            P1=5.D-4
            P2=20.D0
         ENDIF
      ENDIF   
C
      CALL WNTV3(P1,P2,NUTV2N,V,T,EPS,X,IX)
C
      PTV2N=X
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUTV2N(P,T,V)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      NUTV2N=VPT2N(P,T)-V
C
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION PTV5N(T,V)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision NUTV5N
      EXTERNAL NUTV5N
C
      EPS=1.D-6
C
      VZ1=VPT5N(3.D0,T)
      VZ2=VPT5N(0.1D0,T)
      IF (V .LE. VZ1) THEN
         P1=3.D0
         P2=10.D0
      ELSE IF (V .LE. VZ2) THEN
         P1=0.1D0
         P2=3.D0
      ELSE
         P1=5.D-4
         P2=0.1D0
      ENDIF
C
      CALL WNTV3(P1,P2,NUTV5N,V,T,EPS,X,IX)
C
      PTV5N=X
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUTV5N(P,T,V)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      NUTV5N=VPT5N(P,T)-V
C
      END
C
C*******************************************************************
      SUBROUTINE WNTV3(XA,XB,F,P,T,EPS,X,IX)
C*******************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      EXTERNAL F
      X1=XA
      F1=F(X1,T,P)
      X3=XB
      F3=F(X3,T,P)
      IX=0
C++ SCHLEIFE ZUR ERMITTLUNG DER NULLSTELLE +++++++++++++++++++++++++++
      DO 100 I=1,40
      IF( F1 .NE. F3 ) THEN
       X=X1+(X3-X1)*F1/(F1-F3)
      ELSE
       IX=3
       GOTO 999
      END IF
      IF(X .LT. 0.D0) X=(X1+X3)/2.D0
      IF(DABS(X) .LT. 1.D-8) THEN
        IF(DABS(X-X1).LT.EPS) RETURN
      ELSE
        IF(DABS((X-X1)/X).LT.EPS) RETURN
      END IF
      F2=F(X,T,P)
      X2=X1-(X1-X3)/2.D0
      P1=F2*F1
      P3=F2*F3
      IF ((P1 .LT. 0.D0) .AND. (P3 .LT. 0.D0)) THEN
        P1=DABS(P1)
        P3=DABS(P3)
      ENDIF
      IF(P1 .LE. P3) THEN
        X3=X1
        F3=F1
      END IF
      X1=X
      F1=F2
      IF((X2-X3)*(X2-X1) .GE. 0.D0) GOTO 100
c      X=(X1+X3)/2.D0
      IF ((DABS(F1/F3)) .LE. (0.3D0)) THEN
         X=X1+DABS(F1/F3)*1.5D0*(X3-X1)
         IREM=1
      ELSE IF ((DABS(F1/F3)) .GE. (3.D0)) THEN
         IREM=2
      ELSE
         X=(X1+X3)/2.D0
         IREM=0
      ENDIF
      F2=F(X,T,P)
      IF ((((F2*F1) .GE. 0.D0) .AND. (IREM .EQ. 1)) .OR.
     1    (((F2*F3) .GE. 0.D0) .AND. (IREM .EQ. 2))) THEN
         X=(X1+X3)/2.D0
         F2=F(X,T,P)
      ENDIF
      IF(F2*F1 .LE. F2*F3) THEN
        X3=X1
        F3=F1
      END IF
      X1=X
      F1=F2
 100  CONTINUE
C++ ENDE DER SCHLEIFE ++++++++++++++++++++++++++++++++++++++++++++++++
      IX=1
 999  RETURN
      END
C
C***********************************************************************
      SUBROUTINE REGSHV(H,V,IREG)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      COMMON/CCRIT/VC,SC,HC
      COMMON/CGR/PGR,TGR13
C
      double precision HOLD, VOLD
      SAVE HOLD, VOLD
      INTEGER IREGOLD
      SAVE IREGOLD
      data HOLD, VOLD, IREGOLD / -1.D-6, -1.D+0, -1 /
C
      IF ((DABS(H-HOLD).LT. 1.D-6) .AND.
     *   (DABS(V-VOLD).LT. 1.D-6)) THEN
         IREG=IREGOLD
         GOTO 999
      END IF
C
      IF ((V .LT. 9.566869D-4) .OR. (H .LT. -0.05D0) .OR.
     &    (H .GE. 7376.98030619238D0) .OR. 
     &    (V .GE. 2098.23611660241d0)) THEN
         IREG=0
         GOTO 999
      ENDIF
C       
        IREG=0
C
      IF ((H .LT. HC) .AND. (V .LE. VC)) THEN
         PS=PSATH1N(H)
         TS=TSATPN(PS)
         IF (TS .LE. 623.15D0) THEN
            V1=VPT1N(PS,TS)
         ELSE
            CALL FSATP(DV,DL,TOUT,TS,PS)
            V1=1.D0/DL
         ENDIF   
         IF (V .GT. V1) THEN
            IREG=9
         ELSEIF (H .LE. 1553.922503371486D0) THEN
            TN=TPH1N(100.D0,H)
            VN=VPT1N(102.D0,TN)
            IF (V .GE. VN) THEN
               IREG=1
            ELSE
               IREG=0
            ENDIF
         ELSEIF (H .GE. 1670.858218274385D0) THEN
            TN=TPH3N(100.D0,H)
            VN=VPT3N(102.D0,TN)
            IF (V .GE. VN) THEN
               IREG=3
            ELSE
               IREG=0
            ENDIF
         ELSE
C
ccc        Grenze 1/3
C      
            HG1=((((-1.2873623355D18*V+1.1128250625D16)*V-
     &          3.87764700168D13)*V+6.7986880971D10)*V-5.956307233D7)*V+
     &          22275.4105497D0
C
            IF (H .LT. HG1) THEN
               TN=TPH1N(100.D0,H)
               VN=VPT1N(102.D0,TN)
               IF (V .GE. VN) THEN
                  IREG=1
               ELSE
                  IREG=0
               ENDIF
            ELSE
C      
               HG3=(((((-6.6365860607D20*V+4.370326894D18)*V-
     &             8.808032455D15)*V-16.43373971D11)*V+2.947556851D10)*V
     &             -3.851044435D7)*V+17547.8715D0
C
               PHF=PHV3N(H,V)
               IF (H .GT. HG3) THEN
                  IF (PHF .LE. 102.D0) THEN
                     IREG=3
                  ELSE
                     IREG=0
                  ENDIF
               ELSE
                  THF=TPV3N(PHF,V)
                  IF ((THF .LE. 623.15D0) .AND. (PHF .LE. 102.D0)) THEN
                     IREG=1
                  ELSEIF (PHF .LE. 102.D0) THEN
                     IREG=3
                  ELSE
                     IREG=0
                    ENDIF
               ENDIF
            ENDIF  
ccc
         ENDIF
      ELSEIF ((H .GE. HC) .AND. (V .LE. VC)) THEN
         IF (H .LE. 2671.558359836704D0) THEN
            TN=TPH3N(100.D0,H)
            VN=VPT3N(102.D0,TN)
            IF (V .GE. VN) THEN
               IREG=3
            ELSE   
               IREG=0
            ENDIF
         ELSE
ccc      Grenze 2/3 (zwischen 63.66 und 100)   
C            H100=2812.94206d0
            V100=0.25847185d-2
            H90O=2855.22611d0
            H90U=2688.84198d0
            V90=0.268503154d-2
            H80O=2820.83893d0
            H80U=2643.34656d0
            V80=0.280973554d-2
            H70O=2788.70137d0
            H70U=2596.37369d0
            V70=0.29721606d-2
            T60=FB23P(63.66D0)
            H60O=hpt2n(60.d0,T60+15.D0)
            H60U=hpt(60.d0,T60-15.D0,3)
            V60=vpt2n(63.66d0,T60)
            IF (V .LT. V100) THEN
               IREG=0
            ELSEIF ((V .GE. V100) .AND. (V .LE. V90)) THEN
               IF (H .GE. H90O) THEN
                  TG=TPH2N(100.D0,H)
                  VG=VPT2N(102.D0,TG)
                  IF (V .GT. VG) THEN
                     IREG=2
                  ELSE
                     IREG=0
                  ENDIF
               ELSEIF (H .LE. H90U) THEN
                  IREG=3
               ELSE
                  PN=PHV3N(H,V)
                  TN=TPH3N(PN,H)
                  TB=FB23P(PN)
                  IF(TN .LE. TB) THEN
                     IREG=3
                  ELSE
                     TG=TPH2N(100.D0,H)
                     VG=VPT2N(102.D0,TG)
                     IF (V .GT. VG) THEN
                        IREG=2
                     ELSE
                        IREG=0
                     ENDIF
                  ENDIF
               ENDIF
            ELSEIF ((V .GT. V90) .AND. (V .LE. V80)) THEN
               IF (H .GE. H80O) THEN
                  TG=TPH2N(100.D0,H)
                  VG=VPT2N(102.D0,TG)
                  IF (V .GT. VG) THEN
                     IREG=2
                  ELSE
                     IREG=0
                  ENDIF
               ELSEIF (H .LE. H80U) THEN
                  IREG=3
               ELSE
                  PN=PHV3N(H,V)
                  TN=TPH3N(PN,H)
                  TB=FB23P(PN)
                  IF(TN .LE. TB) THEN
                     IREG=3
                  ELSE
                     TG=TPH2N(100.D0,H)
                     VG=VPT2N(102.D0,TG)
                     IF (V .GT. VG) THEN
                        IREG=2
                     ELSE
                        IREG=0
                     ENDIF
                  ENDIF
               ENDIF
            ELSEIF ((V .GT. V80) .AND. (V .LE. V70)) THEN
               IF (H .GE. H70O) THEN
                  TG=TPH2N(100.D0,H)
                  VG=VPT2N(102.D0,TG)
                  IF (V .GT. VG) THEN
                     IREG=2
                  ELSE
                     IREG=0
                  ENDIF
               ELSEIF (H .LE. H70U) THEN
                  IREG=3
               ELSE
                  PN=PHV3N(H,V)
                  TN=TPH3N(PN,H)
                  TB=FB23P(PN)
                  IF(TN .LE. TB) THEN
                     IREG=3
                  ELSE
                     TG=TPH2N(100.D0,H)
                     VG=VPT2N(102.D0,TG)
                     IF (V .GT. VG) THEN
                        IREG=2
                     ELSE
                        IREG=0
                     ENDIF
                  ENDIF
               ENDIF
            ELSEIF ((V .GT. V70) .AND. (V .LE. V60)) THEN
               IF (H .GE. H60O) THEN
                  TG=TPH2N(100.D0,H)
                  VG=VPT2N(102.D0,TG)
                  IF (V .GT. VG) THEN
                     IREG=2
                  ELSE
                     IREG=0
                  ENDIF
               ELSEIF (H .LE. H60U) THEN
                  IREG=3
               ELSE
                  PN=PHV3N(H,V)
                  TN=TPH3N(PN,H)
                  TB=FB23P(PN)
                  IF(TN .LE. TB) THEN
                     IREG=3
                  ELSE
                     TG=TPH2N(100.D0,H)
                     VG=VPT2N(102.D0,TG)
                     IF (V .GT. VG) THEN
                        IREG=2
                     ELSE
                        IREG=0
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
cccc
         ENDIF
      ELSE
         IF (V .LE. 205.997459488978d0) THEN
              PS=PSATV2N(V)
            TS=TSATPN(PS)
            IF (TS .LE. 623.15D0) THEN
               H2=HPT2N(PS,TS)
            ELSE
               CALL FSATP(DV,DL,TOUT,TS,PS)
               V2=1.D0/DV
               H2=HVT3N(V2,TS)
            ENDIF
           ELSE
              H2=HC
           ENDIF
C
         IF (H .LE. H2) THEN
            IREG=9
         ELSEIF ((H .LE. 2671.558359836704D0) .AND. 
     &    (V .LE. 8.800931931575926D-3)) THEN
ccc      Grenze 2/3 (zwischen 63.66 und Sat)
            T70=FB23P(63.66D0)
            H70O=hpt2n(70.d0,T70+15.D0)
            H70U=hpt(70.d0,T70-15.D0,3)
            V70=vpt2n(63.66d0,T70)
            H60O=2760.81633d0
            H60U=2547.83462d0
            V60=0.319786479d-2
            H50O=2740.88557d0
            H50U=2497.87983d0
            V50=0.353809389d-2
            H40O=2735.18340d0
            H40U=2447.30113d0
            V40=0.412070496d-2
            H30O=2748.86452d0
            H30U=2395.83710d0
            V30=0.530032498d-2
C            H20O=2759.52306d0
C            H20U=1759.88222d0
C            V20=0.787818006d-2
            VSAT=vpt2n(PGR,623.15d0)
            IF ((V .GE. V70) .AND. (V .LE. V60)) THEN
               IF (H .GE. H60O) THEN
                  IREG=2
               ELSEIF (H .LE. H60U) THEN
                  IREG=3
               ELSE
                  PN=PHV3N(H,V)
                  TN=TPH3N(PN,H)
                  TB=FB23P(PN)
                  IF(TN .LE. TB) THEN
                     IREG=3
                  ELSE
                     IREG=2
                  ENDIF
               ENDIF
            ELSEIF ((V .GT. V60) .AND. (V .LE. V50)) THEN
               IF (H .GE. H50O) THEN
                  IREG=2
               ELSEIF (H .LE. H50U) THEN
                  IREG=3
               ELSE
                  PN=PHV3N(H,V)
                  TN=TPH3N(PN,H)
                  TB=FB23P(PN)
                  IF(TN .LE. TB) THEN
                     IREG=3
                  ELSE
                     IREG=2
                  ENDIF
               ENDIF
            ELSEIF ((V .GT. V50) .AND. (V .LE. V40)) THEN
               IF (H .GE. H40O) THEN
                  IREG=2
               ELSEIF (H .LE. H40U) THEN
                  IREG=3
               ELSE
                  PN=PHV3N(H,V)
                  TN=TPH3N(PN,H)
                  TB=FB23P(PN)
                  IF(TN .LE. TB) THEN
                     IREG=3
                  ELSE
                     IREG=2
                  ENDIF
               ENDIF
            ELSEIF ((V .GT. V40) .AND. (V .LE. V30)) THEN
               IF (H .GE. H30O) THEN
                  IREG=2
               ELSEIF (H .LE. H30U) THEN
                  IREG=3
               ELSE
                  PN=PHV3N(H,V)
                  TN=TPH3N(PN,H)
                  TB=FB23P(PN)
                  IF(TN .LE. TB) THEN
                     IREG=3
                  ELSE
                     IREG=2
                  ENDIF
               ENDIF
            ELSEIF ((V .GT. V30) .AND. (V .LE. VSAT)) THEN
               HB=HBV(V)
               IF (H .GT. HB) THEN
                  IREG=2
               ELSE
                  IREG=3
               ENDIF
            ENDIF
         ELSEIF (H .LE. 3715.188943537408D0) THEN
            TN=TPH2C(100.D0,H)
            VN=VPT2N(102.D0,TN)
            IF (V .GE. VN) THEN
               IREG=2
            ELSE   
               IREG=0
            ENDIF
         ELSEIF (V .LE. 0.433550765d-2) THEN
              IREG=0
         ELSEIF (V .LE. 990.572347145D0) THEN
            D=1.D0/V
              VN=D+2500.D0
C
            A5 = 2.99193993954271D-11
            A4 = -4.01914437088434D-07
            A3 = 2.15845404507264D-03
            A2 = -5.79116740809968D0
            A1 = 7.75859470381492D3
            A0 = -4.14540785822924D6
C
            B6 = -6.94169012547241D-14
            B5 = 8.13115207549399D-11
            B4 = -3.1227506664461D-08
            B3 = 5.20112279349633D-06
            B2 = -3.81153980483174D-04
            B1 = 1.01231272408059D-02
            B0 = 1.3869135097795d-3
C
            HGRENZ = ((((A5*VN+A4)*VN+A3)*VN+A2)*VN+A1)*VN+A0
     &       + (((((B6*D+B5)*D+B4)*D+B3)*D+B2)*D+B1)*D+B0
C
            IF (H .LT. HGRENZ) THEN
                 TN=TPH2N(5.D-4,H)
                 VN=VPT2N(4.99D-4,TN)
                 IF (V .LE. VN) THEN
                            IREG=2
                 ELSE
                    IREG=0
                 ENDIF
              ELSEIF (V .LT. 4.862419121720567D-2) THEN
                 IREG=0
              ELSEIF (H .LT. 4160.67903103581D0) THEN
               TN=TPH5N(10.D0,H)
                 VN=VPT5N(10.2D0,TN)
               IF (V .GE. VN) THEN
                  IREG=5
               ELSE
                  IREG=0
               ENDIF
              ELSEIF (H .LT. 7374.75171191195D0) THEN
                 TNO=TPH5N(10.D0,H)
                 VNO=VPT5N(10.D0,TNO)
                 TNU=TPH5N(5.D-4,H)
                 VNU=VPT5N(4.99D-4,TNU)
               IF ((V .GE. VNO) .AND. (V .LE. VNU)) THEN
                  IREG=5
               ELSE
                  IREG=0
               ENDIF
              ELSEIF (H .LT. 7376.98030619238D0) THEN
                 PN=PHV5N(H,V)
                 TN=TPH5N(PN,H)
               IF (TN .LE. 2273.2D0) THEN
                  IREG=5
               ELSE
                  IREG=0
               ENDIF
              ELSE
                 IREG=0
              ENDIF
           ELSE
c            PN=PTV5N(2273.15d0,V)
c             HN=HPT5N(PN,2273.15D0)
            IF (H .GT. 7376.98030619238D0) THEN
               IREG=0
              ELSE
               TN=TPH5N(5.D-4,H)
               VN=VPT5N(4.99D-4,TN)
                 IF (V .LE. VN) THEN
                    IREG=5
                 ELSE
                    IREG=0
                 ENDIF
            ENDIF
           ENDIF      
      ENDIF
C
999   CONTINUE
C
      HOLD=H
      VOLD=V
      IREGOLD=IREG
      RETURN   
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION PHV(H,V,IREG)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      IF (IREG .EQ. 0) THEN
         CALL REGSHV(H,V,IREG)
      ENDIF
C  
      IF (IREG .EQ. 1) THEN
         PHV=PHV1N(H,V)
      ELSE IF (IREG .EQ. 2) THEN
         PHV=PHV2N(H,V)
      ELSE IF (IREG .EQ. 3) THEN
         PHV=PHV3N(H,V)
      ELSE IF (IREG .EQ. 5) THEN
         PHV=PHV5N(H,V)
      ELSE IF (IREG .EQ. 9) THEN
         PHV=PHV9N(H,V)
      ENDIF
C
      END
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION TPH1C(P,H)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      T1=TPH1N(P,H)
      H1=HPT1N(P,T1)
      HN=2.D0*H-H1
      cp=cppt1n(p,t1)
      T2=T1+(h-h1)/cp
      H2=HPT1N(P,T2)
      TPH1C=T2+(T1-T2)*(H-H2)/(H1-H2)
C
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION PHV1N(H,V)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision NUHV1N
      EXTERNAL NUHV1N
C
      EPS=1.D-6
C
      TZ1=TPH1N(60.D0,H)
      VZ1=VPT1N(60.D0,TZ1)
      TZ2=TPH1N(25.D0,H)
      VZ2=VPT1N(25.D0,TZ2)
      IF (V .LE. VZ1) THEN
         P1=60.D0
         P2=100.D0
      ELSE IF (V .LE. VZ2) THEN
         P1=25.D0
         P2=60.D0
      ELSE
         P1=PSATTN(TZ2)
         P2=25.D0
      ENDIF
C
      CALL WNVH3(P1,P2,NUHV1N,V,H,EPS,X,IX)
C
      PHV1N=X
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUHV1N(P,H,V)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      T=TPH1C(P,H)
      NUHV1N=VPT1N(P,T)-V
C
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION PSATV2N(V)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision NUV2SATN
      EXTERNAL NUV2SATN
C
      EPS=1.D-9
C
      IF (V .LE. 0.008828D0) THEN
         P1=18.1d0
         P2=22.064d0
      ELSEIF (V .LE. 0.2D0) THEN
         P1=1.d0
         P2=10.d0
      ELSEIF (V .LE. 1.7D0) THEN
         P1=0.1d0
         P2=1.d0
      ELSEIF (V .LE. 15.3D0) THEN
         P1=0.01d0
         P2=0.1d0
      ELSE
         P1=0.0005D0
         P2=0.01D0
      ENDIF
      DUM=1.D0
C
      CALL WNHSAT3(P1,P2,NUV2SATN,DUM,V,EPS,X,IX)
C
      PSATV2N=X
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUV2SATN(P,V,DUM)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      Z=DUM
      T=TSATPN(P)
C
      IF (T .LE. 623.15D0) THEN
         VN=VPT2N(P,T)
      ELSE
         CALL FSATP(DV,DL,TOUT,T,P)
         VN=1.D0/DV
      ENDIF
C
      NUV2SATN=VN-V
C
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION HBV(V)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision NUVBN
      EXTERNAL NUVBN
C
      EPS=1.D-6
C
      DUM = -1.D+0
      P1=PSATTN(623.15D0)
      P2=20.D0
C
      CALL WNHSAT3(P1,P2,NUVBN,DUM,V,EPS,X,IX)
C
      TB=FB23P(x)
      HBV=HPT2N(X,TB)
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUVBN(P,V,DUM)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      Z=DUM
      PG=PSATTN(623.15D0)
      IF (P .LT. PG) THEN
         PN=PG
      ELSE   
         PN=P
      ENDIF
      T=FB23P(PN)
      NUVBN=VPT2N(PN,T)-V
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION PHV2N(H,V)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision NUHV2N
      EXTERNAL NUHV2N
C
      EPS=1.D-6
C
      IF (H .GT. 2850D0) THEN
         IF (V .LT. 1.D-2) THEN
            P1=40.D0
            P2=70.D0
         ELSEIF (V .LT. 2.D-2) THEN
            P1=10.D0
            P2=20.D0
         ELSEIF (V .LT. 2.D-1) THEN
            P1=1.D0
            P2=10.D0
         ELSEIF (V .LT. 2.D0) THEN
            P1=1.D-1
            P2=1.D0
         ELSEIF (V .LT. 2.D1) THEN
            P1=1.D-2
            P2=1.D-1
         ELSEIF (V .LT. 2.D2) THEN
            P1=1.D-3
            P2=1.D-2
         ELSE
            P1=5.D-4
            P2=1.D-3
         ENDIF
      ELSEIF (V .LT. 1.8D-2) THEN
         P1=10.D0
         P2=15.D0
      ELSEIF (V .LT. 1.95D-1) THEN
         P1=1.D0
         P2=10.D0
      ELSEIF (V .LT. 1.7D0) THEN
         P1=1.D-1
         P2=1.D0
      ELSEIF (V .LT. 1.467D1) THEN
         P1=1.D-2
         P2=1.D-1
      ELSEIF (V .LT. 1.292D2) THEN
         P1=1.D-3
         P2=1.D-2
      ELSE
         P1=5.D-4
         P2=1.D-3
      ENDIF
C
      CALL WNVH3(P1,P2,NUHV2N,V,H,EPS,X,IX)
C
      PHV2N=X
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUHV2N(P,H,V)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      T=TPH2N(P,H)
      NUHV2N=VPT2N(P,T)-V
C
      END
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION TPH2C(P,H)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      T1=TPH2N(P,H)
      H1=HPT2N(P,T1)
      CP=CPPT2N(P,T1)
      T2=T1+(h-h1)/cp
      H2=HPT2N(P,T2)
      TPH2C=T2+(T1-T2)*(H-H2)/(H1-H2)
C
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION PHV3N(H,V)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision NUHV3N
      EXTERNAL NUHV3N
C
      EPS=1.D-6
C
      T1=650.d0
      T2=800.D0
C
      CALL WNVH3(T1,T2,NUHV3N,H,V,EPS,X,IX)
C
      PHV3N=PVT3N(V,X)
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUHV3N(T,V,H)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      NUHV3N=HVT3N(V,T)-H
C
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION PHV5N(H,V)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision NUHV5N
      EXTERNAL NUHV5N
C
      EPS=1.D-6
C
      IF (V .LE. 5.D0) THEN
           P1=0.1D0
           P2=10.D0
        ELSE
           P1=0.001d0
           P2=0.21d0
        ENDIF
C
      CALL WNVH3(P1,P2,NUHV5N,V,H,EPS,X,IX)
C
      PHV5N=X
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUHV5N(P,H,V)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      T=TPH(P,H,5)
      NUHV5N=VPT5N(P,T)-V
C
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION PHV9N(H,V)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision NUHV9N
      EXTERNAL NUHV9N
      COMMON/CCRIT/VC,SC,HC
C
      EPS=1.D-6
C
      IF (H .LE. HC) THEN
         PG=PSATH1N(H)
      ELSE
         PG=PSATV2N(V)
      ENDIF
      IF (PG .GE. 1.D-3) THEN
           P1=1.D-3
           P2=PG
        ELSE
           P1=6.11657D-4
           P2=PG
        ENDIF
C
      CALL WNVH3(P1,P2,NUHV9N,V,H,EPS,X,IX)
C
      PHV9N=X
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUHV9N(P,H,V)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      IF (P .GT. 22.064D0) THEN
         PN=22.06D0
      ELSEIF (P .LT. 5.D-4) THEN
         PN=5.D-4
      ELSE
         PN=P
      ENDIF
      T=TSATPN(PN)
      IF (T .LE. 623.15D0) THEN
         H1=HPT1N(PN,T)
         H2=HPT2N(PN,T)
         V1=VPT1N(PN,T)
         V2=VPT2N(PN,T)
      ELSE
         CALL FSATP(DV,DL,TOUT,T,PN)
         V1=1.D0/DL
         V2=1.D0/DV
         H1=HVT3N(V1,T)
         H2=HVT3N(V2,T)
      ENDIF
      NUHV9N=H1+(V-V1)/(V2-V1)*(H2-H1)-H
C
      END
C
C*******************************************************************
      SUBROUTINE WNVH3(XA,XB,F,P,T,EPS,X,IX)
C*******************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      EXTERNAL F
      X1=XA
      F1=F(X1,T,P)
      X3=XB
      F3=F(X3,T,P)
      IX=0
C++ SCHLEIFE ZUR ERMITTLUNG DER NULLSTELLE +++++++++++++++++++++++++++
      DO 100 I=1,40
      IF( F1 .NE. F3 ) THEN
       X=X1+(X3-X1)*F1/(F1-F3)
      ELSE
       IX=3
       GOTO 999
      END IF
      IF(X .LT. 0.D0) X=(X1+X3)/2.D0
      IF(DABS(X) .LT. 1.D-8) THEN
        IF(DABS(X-X1).LT.EPS) RETURN
      ELSE
        IF(DABS((X-X1)/X).LT.EPS) RETURN
      END IF
      F2=F(X,T,P)
      X2=X1-(X1-X3)/2.D0
      P1=F2*F1
      P3=F2*F3
      IF ((P1 .LT. 0.D0) .AND. (P3 .LT. 0.D0)) THEN
        P1=DABS(P1)
        P3=DABS(P3)
      ENDIF
      IF(P1 .LE. P3) THEN
        X3=X1
        F3=F1
      END IF
      X1=X
      F1=F2
      IF((X2-X3)*(X2-X1) .GE. 0.D0) GOTO 100
c      X=(X1+X3)/2.D0
      IF ((DABS(F1/F3)) .LE. (0.3D0)) THEN
         X=X1+DABS(F1/F3)*1.5D0*(X3-X1)
         IREM=1
      ELSE IF ((DABS(F1/F3)) .GE. (3.D0)) THEN
         IREM=2
      ELSE
         X=(X1+X3)/2.D0
         IREM=0
      ENDIF
      F2=F(X,T,P)
      IF ((((F2*F1) .GE. 0.D0) .AND. (IREM .EQ. 1)) .OR.
     1    (((F2*F3) .GE. 0.D0) .AND. (IREM .EQ. 2))) THEN
         X=(X1+X3)/2.D0
         F2=F(X,T,P)
      ENDIF
      IF(F2*F1 .LE. F2*F3) THEN
        X3=X1
        F3=F1
      END IF
      X1=X
      F1=F2
 100  CONTINUE
C++ ENDE DER SCHLEIFE ++++++++++++++++++++++++++++++++++++++++++++++++
      IX=1
 999  RETURN
      END
C
C*******************************************************************
      SUBROUTINE WNHSAT3(XA,XB,F,P,T,EPS,X,IX)
C*******************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      EXTERNAL F
      X1=XA
      F1=F(X1,T,P)
      X3=XB
      F3=F(X3,T,P)
      IX=0
C++ SCHLEIFE ZUR ERMITTLUNG DER NULLSTELLE +++++++++++++++++++++++++++
      DO 100 I=1,40
      IF( F1 .NE. F3 ) THEN
       X=X1+(X3-X1)*F1/(F1-F3)
      ELSE
       IX=3
       GOTO 999
      END IF
      IF(X .LT. 0.D0) X=(X1+X3)/2.D0
      IF(DABS(X) .LT. 1.D-8) THEN
        IF(DABS(X-X1).LT.EPS) RETURN
      ELSE
        IF(DABS((X-X1)/X).LT.EPS) RETURN
      END IF
      F2=F(X,T,P)
      X2=X1-(X1-X3)/2.D0
      P1=F2*F1
      P3=F2*F3
      IF ((P1 .LT. 0.D0) .AND. (P3 .LT. 0.D0)) THEN
        P1=DABS(P1)
        P3=DABS(P3)
      ENDIF
      IF(P1 .LE. P3) THEN
        X3=X1
        F3=F1
      END IF
      X1=X
      F1=F2
      IF((X2-X3)*(X2-X1) .GE. 0.D0) GOTO 100
c      X=(X1+X3)/2.D0
      IF ((DABS(F1/F3)) .LE. (0.3D0)) THEN
         X=X1+DABS(F1/F3)*1.5D0*(X3-X1)
         IREM=1
      ELSE IF ((DABS(F1/F3)) .GE. (3.D0)) THEN
         IREM=2
      ELSE
         X=(X1+X3)/2.D0
         IREM=0
      ENDIF
      F2=F(X,T,P)
      IF ((((F2*F1) .GE. 0.D0) .AND. (IREM .EQ. 1)) .OR.
     1    (((F2*F3) .GE. 0.D0) .AND. (IREM .EQ. 2))) THEN
         X=(X1+X3)/2.D0
         F2=F(X,T,P)
      ENDIF
      IF(F2*F1 .LE. F2*F3) THEN
        X3=X1
        F3=F1
      END IF
      X1=X
      F1=F2
 100  CONTINUE
C++ ENDE DER SCHLEIFE ++++++++++++++++++++++++++++++++++++++++++++++++
      IX=1
 999  RETURN
      END
C
C***********************************************************************
      SUBROUTINE REGSSV(S,V,IREG)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      COMMON/CCRIT/VC,SC,HC
      COMMON/CGR/PGR,TGR13
C
      double precision SOLD, VOLD
      SAVE SOLD, VOLD
      INTEGER IREGOLD
      SAVE IREGOLD
      data SOLD, VOLD, IREGOLD / -1.D-6, -1.D+0, -1 /
C
      IF ((DABS(S-SOLD).LT. 1.D-6) .AND.
     *   (DABS(V-VOLD).LT. 1.D-6)) THEN
         IREG=IREGOLD
         GOTO 999
      END IF
C
      IF ((V .LT. 9.566869D-4) .OR. (S .LT. -0.5D0) .OR.
     &    (S .GE. 13.99764756652987D0) .OR. (V .GE. 2099)) THEN
         IREG=0
         GOTO 999
      ENDIF
C
        IREG=0
C
      IF ((S .LT. SC) .AND. (V .LE. VC)) THEN
         PS=PSATS1N(S)
         TS=TSATPN(PS)
         IF (TS .LE. 623.15D0) THEN
            V1=VPT1N(PS,TS)
         ELSE
            CALL FSATP(DV,DL,TOUT,TS,PS)
            V1=1.D0/DL
         ENDIF   
         IF (V .GT. V1) THEN
            IREG=9
ccc
         ELSEIF (S .LE. 3.397782954701892D0) THEN
            TN=TPS1N(100.D0,S)
            VN=VPT1N(102.D0,TN)
            IF (V .GE. VN) THEN
               IREG=1
            ELSE
               IREG=0
            ENDIF
ccc
         ELSEIF (S .GE. 3.778281339544279D0) THEN
            TN=TPS3N(100.D0,S)
            VN=VPT3N(102.D0,TN)
            IF (V .GE. VN) THEN
               IREG=3
            ELSE
               IREG=0
            ENDIF
         ELSE
ccc        Grenze 1/3
            VN=(V*1.D3)-1.3D0
            SG1=-1.74650514276785d-2+vn*(1.25791078740034d0+vn*(
     &       -1.02344008564254d0+vn*4.68814970415902d-1))+3.4d0      
C
            IF (S .LT. SG1) THEN
               TN=TPS1N(100.D0,S)
               VN=VPT1N(102.D0,TN)
               IF (V .GE. VN) THEN
                  IREG=1
               ELSE
                  IREG=0
               ENDIF
            ELSE
               SG3=SG1+5.d-3
               PHF=PSV3N(S,V)
               IF (S .GT. SG3) THEN
                            IF (PHF .LE. 102.D0) THEN
                     IREG=3
                  ELSE
                     IREG=0
                  ENDIF
               ELSE
                  THF=TPS3N(PHF,S)
                  IF ((THF .LE. 623.15D0) .AND. (PHF .LE. 102.D0)) THEN
                     IREG=1
                  ELSEIF (PHF .LE. 102.D0) THEN
                     IREG=3
                  ELSE
                     IREG=0
                    ENDIF
               ENDIF
            ENDIF  
         ENDIF
C
      ELSEIF ((S .GE. SC) .AND. (V .LE. VC)) THEN
         IF (S .LE. 5.D0) THEN
            TN=TPS3N(100.D0,S)
            VN=VPT3N(102.D0,TN)
            IF (V .GE. VN) THEN
               IREG=3
            ELSE   
               IREG=0
            ENDIF
         ELSE
            IF (S .GT. 5.435D0) THEN
               IREG=0
               GOTO 999
            ENDIF
ccc      Grenze 2/3 (zwischen 63.66 und 100)   
C            S100=5.09796903d0
            V100=0.25847185d-2
C            S90=5.08243607d0
            S90O=5.17841746d0
            S90U=4.98163547d0
            V90=0.268503154d-2
C            S80=5.06786979d0
            S80O=5.17841746d0
            S80U=4.98163647d0
            V80=0.250973554d-2
C            S70=5.05590949d0
            S70O=5.17088967d0
            S70U=4.93249626d0
            V70=0.297216060d-2
            T60=FB23P(63.66D0)
            S60O=spt2n(60.d0,T60+15.D0)
            S60U=spt(60.d0,T60-15.D0,3)
            V60=vpt2n(63.66d0,T60)
              IF (V .LT. 0.00236D0) THEN
               IREG=0
              ELSEIF (V .LT. V100) THEN
               TG=TPS3N(100.d0,S)
                 VG=VPT3N(100.d0,TG)
                         IF (V .LT. VG) THEN
                            IREG=0
                 ELSE
                    IREG=3
                 ENDIF
            ELSEIF ((V .GE. V100) .AND. (V .LE. V90)) THEN
               IF (S .GE. S90O) THEN
                  TG=TPS2N(100.D0,S)
                  VG=VPT2N(102.D0,TG)
                  IF (V .GT. VG) THEN
                     IREG=2
                  ELSE
                     IREG=0
                  ENDIF
               ELSEIF (S .LE. S90U) THEN
                  IREG=3
               ELSE
                  PN=PSV3N(S,V)
                  TN=TPS3N(PN,S)
                  TB=FB23P(PN)
                  IF(TN .LE. TB) THEN
                     IREG=3
                  ELSE
                     TG=TPS2N(100.D0,S)
                     VG=VPT2N(102.D0,TG)
                     IF (V .GT. VG) THEN
                        IREG=2
                     ELSE
                        IREG=0
                     ENDIF
                  ENDIF
               ENDIF
            ELSEIF ((V .GT. V90) .AND. (V .LE. V80)) THEN
               IF (S .GE. S80O) THEN
                  TG=TPS2N(100.D0,S)
                  VG=VPT2N(102.D0,TG)
                  IF (V .GT. VG) THEN
                     IREG=2
                  ELSE
                     IREG=0
                  ENDIF
               ELSEIF (S .LE. S80U) THEN
                  IREG=3
               ELSE
                  PN=PSV3N(S,V)
                  TN=TPS3N(PN,S)
                  TB=FB23P(PN)
                  IF(TN .LE. TB) THEN
                     IREG=3
                  ELSE
                     TG=TPS2N(100.D0,S)
                     VG=VPT2N(102.D0,TG)
                     IF (V .GT. VG) THEN
                        IREG=2
                     ELSE
                        IREG=0
                     ENDIF
                  ENDIF
               ENDIF
            ELSEIF ((V .GT. V80) .AND. (V .LE. V70)) THEN
               IF (S .GE. S70O) THEN
                  TG=TPS2N(100.D0,S)
                  VG=VPT2N(102.D0,TG)
                  IF (V .GT. VG) THEN
                     IREG=2
                  ELSE
                     IREG=0
                  ENDIF
               ELSEIF (S .LE. S70U) THEN
                  IREG=3
               ELSE
                  PN=PSV3N(S,V)
                  TN=TPS3N(PN,S)
                  TB=FB23P(PN)
                  IF(TN .LE. TB) THEN
                     IREG=3
                  ELSE
                     TG=TPS2N(100.D0,S)
                     VG=VPT2N(102.D0,TG)
                     IF (V .GT. VG) THEN
                        IREG=2
                     ELSE
                        IREG=0
                     ENDIF
                  ENDIF
               ENDIF
            ELSEIF ((V .GT. V70) .AND. (V .LE. V60)) THEN
               IF (S .GE. S60O) THEN
                  TG=TPS2N(100.D0,S)
                  VG=VPT2N(102.D0,TG)
                  IF (V .GT. VG) THEN
                     IREG=2
                  ELSE
                     IREG=0
                  ENDIF
               ELSEIF (S .LE. S60U) THEN
                  IREG=3
               ELSE
                  PN=PSV3N(S,V)
                  TN=TPS3N(PN,S)
                  TB=FB23P(PN)
                  IF(TN .LE. TB) THEN
                     IREG=3
                  ELSE
                     TG=TPS2N(100.D0,S)
                     VG=VPT2N(102.D0,TG)
                     IF (V .GT. VG) THEN
                        IREG=2
                     ELSE
                        IREG=0
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
cccc
         ENDIF
      ELSE
         PS=PSATV2N(V)
         TS=TSATPN(PS)
         IF (TS .LE. 623.15D0) THEN
            S2=SPT2N(PS,TS)
         ELSE
            CALL FSATP(DV,DL,TOUT,TS,PS)
            V2=1.D0/DV
            S2=SVT3N(V2,TS)
         ENDIF
         IF (S .LE. S2) THEN
            IREG=9
         ELSEIF ((S .LE. 5.265D0) .AND. 
     &    (V .LE. 8.800931931575926D-3)) THEN
ccc      Grenze 2/3 (zwischen 63.66 und Sat)
            T70=FB23P(63.66D0)
            S70O=spt2n(70.d0,T70+15.D0)
            S70U=spt(70.d0,T70-15.D0,3)
            V70=vpt2n(63.66d0,T70)
            S60O=5.17778011d0
            S60U=4.90636497d0
            V60=0.319786479d-2
            S50O=5.19965504d0
            S50U=4.87990848d0
            V50=0.353809389d-2
            S40O=5.24885292d0
            S40U=4.85519681d0
            V40=0.412070496d-2
            S30O=5.34159153d0
            S30U=4.83419222d0
            V30=0.530032498d-2
C            S20O=5.46677932d0
C            S20U=3.90986399d0
C            V20=0.787818006d-2
            VSAT=vpt2n(PGR,623.15d0)
            IF ((V .GE. V70) .AND. (V .LE. V60)) THEN
               IF (S .GE. S60O) THEN
                  IREG=2
               ELSEIF (S .LE. S60U) THEN
                  IREG=3
               ELSE
                  PN=PSV3N(S,V)
                  TN=TPS3N(PN,S)
                  TB=FB23P(PN)
                  IF(TN .LE. TB) THEN
                     IREG=3
                  ELSE
                     IREG=2
                  ENDIF
               ENDIF
            ELSEIF ((V .GT. V60) .AND. (V .LE. V50)) THEN
               IF (S .GE. S50O) THEN
                  IREG=2
               ELSEIF (S .LE. S50U) THEN
                  IREG=3
               ELSE
                  PN=PSV3N(S,V)
                  TN=TPS3N(PN,S)
                  TB=FB23P(PN)
                  IF(TN .LE. TB) THEN
                     IREG=3
                  ELSE
                     IREG=2
                  ENDIF
               ENDIF
            ELSEIF ((V .GT. V50) .AND. (V .LE. V40)) THEN
               IF (S .GE. S40O) THEN
                  IREG=2
               ELSEIF (S .LE. S40U) THEN
                  IREG=3
               ELSE
                  PN=PSV3N(S,V)
                  TN=TPS3N(PN,S)
                  TB=FB23P(PN)
                  IF(TN .LE. TB) THEN
                     IREG=3
                  ELSE
                     IREG=2
                  ENDIF
               ENDIF
            ELSEIF ((V .GT. V40) .AND. (V .LE. V30)) THEN
               IF (S .GE. S30O) THEN
                  IREG=2
               ELSEIF (S .LE. S30U) THEN
                  IREG=3
               ELSE
                  PN=PSV3N(S,V)
                  TN=TPS3N(PN,S)
                  TB=FB23P(PN)
                  IF(TN .LE. TB) THEN
                     IREG=3
                  ELSE
                     IREG=2
                  ENDIF
               ENDIF
            ELSEIF ((V .GT. V30) .AND. (V .LE. VSAT)) THEN
               SB=SBV(V)
               IF (S .GT. SB) THEN
                  IREG=2
               ELSE
                  IREG=3
               ENDIF
            ENDIF
         ELSEIF (S .LE. 6.040483671712074D0) THEN
            TN=TPS2N(100.D0,S)
            VN=VPT2N(102.D0,TN)
            IF (V .GE. VN) THEN
               IREG=2
            ELSE   
               IREG=0
            ENDIF
         ELSEIF (V .LE. 0.433550764d-2) THEN
              IREG=0
         ELSEIF (V .LE. 0.486241912d-1) THEN
            IF (V .LE. 0.016d0) THEN      
               SG2= 4.486222853d0+v*(6.809098338d2+v*(
     &              -114902.d0+v*(12253809.23d0+v*(
     &              -777559798.d0+v*(2.684956351d10-v*
     &              3.880379209d11)))))
              ELSEIF (V .LE. 0.032d0) THEN
               SG2= 5.55967561666583d0+v*(1.7395646919d2+v*(
     &              -1.04843599d4+v*(4.239442439d5+v*(
     &              -1.055768713d7+v*(1.469372027d8-v*
     &              8.753213184d8)))))
              ELSE
               SG2= 6.25414133543196d0+v*(5.15789d1+v*(
     &              -1.0132143d3+v*(1.299331d4+v*(
     &              -9.8634d4+v*(4.02332d5-v*
     &              6.77638d5)))))
            ENDIF
            IF (S .LE. SG2) THEN
                 IREG=2
              ELSE
                 IREG=0
              ENDIF
         ELSEIF (V .LE. 991.d0) THEN
            IF (V .LE. 0.15d0) THEN      
               SG5= 6.25429133543196d0+v*(5.15789d1+v*(
     &              -1.0132143d3+v*(1.299331d4+v*(
     &              -9.8634d4+v*(4.02332d5-v*
     &              6.77638d5)))))
              ELSEIF (V .LE. 1.d0) THEN
                 SG5= 7.22426164703249d0+v*(7.2834d0+v*(
     &              -2.18468d1+v*(4.32224d1+v*(
     &              -5.06037d1+v*(3.1749d1-v*
     &              8.1998d0)))))
              ELSEIF (V .LE. 6.d0) THEN
               SG5= 8.08108299999443d0+v*(1.16361848545d0+v*(
     &              -5.64799988375d-1+v*(1.8078824981d-1+v*(
     &              -3.41426d-2+v*(3.443657d-3-v*
     &              1.4255d-4)))))
              ELSEIF (V .LE. 30.d0) THEN
               SG5= 8.8851239448354d0+v*(2.06557d-1+v*(
     &              -1.82718d-2+v*(1.092615d-3+v*(
     &              -3.9418d-5+v*(7.7405d-7-v*
     &              6.336d-9)))))
              ELSEIF (V .LE. 200.d0) THEN
               SG5= 9.69240982779412d0+v*(3.564061d-2+v*(
     &              -5.325019779d-4+v*(5.2722d-6+v*(
     &              -3.095069d-8+v*(9.74722817d-11-v*
     &              1.2645d-13)))))
              ELSE
               SG5= 1.05033d1+v*(6.2018d-3+v*(
     &              -1.6474d-5+v*(2.9576d-8+v*(
     &              -3.20276d-11+v*(1.8874d-14-v*
     &              4.6353d-18)))))
            ENDIF
C
              IF (S .LT. SG5) THEN
                 SG2=SG5-2.5D-3
                 IF (S .LT. SG2) THEN
                    IREG=2
                 ELSE
                    PHF=PSV5N(S,V)
                    THF=TPS5N(PHF,S)
                    IF (THF .LE. 1073.15D0) THEN
                       IREG=2
                    ELSEIF (PHF .LE. 10.2D0) THEN
                       IREG=5
                    ELSE
                       IREG=0
                    ENDIF
                 ENDIF
              ELSEIF (V .LT. 0.105377995D0) THEN
                 TNO=TPV5N(10.2D0,V)
                 SNO=SPT5N(10.D0,TNO)
                     IF (S .LE. SNO) THEN
                    IREG=5
                 ELSE
                    IREG=0
                 ENDIF
            ELSE
                 PN=PSV5N(S,V)
                 TN=TPS5N(PN,S)
                 IF (TN .LE. 2273.2D0) THEN
                  IREG=5
                 ELSE
                    IREG=0
                 ENDIF
              ENDIF
           ELSEIF (S .LE. 13.99764756652987D0) THEN
              PN=PSV5N(S,V)
              TN=TPS5N(PN,S)
              IF ((TN .LE. 2273.2D0) .AND. (PN .GE. 5.D-4)) THEN
               IREG=5
              ELSE
                     IREG=0
              ENDIF
         ELSE
              IREG=0
           ENDIF
      ENDIF
C
999   CONTINUE
C
      SOLD=S
      VOLD=V
      IREGOLD=IREG
C   
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION PSV(S,V,IREG)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      IF (IREG .EQ. 0) THEN
         CALL REGSSV(S,V,IREG)
      ENDIF
C  
      IF (IREG .EQ. 1) THEN
         PSV=PSV1N(S,V)
      ELSE IF (IREG .EQ. 2) THEN
         PSV=PSV2N(S,V)
      ELSE IF (IREG .EQ. 3) THEN
         PSV=PSV3N(S,V)
      ELSE IF (IREG .EQ. 5) THEN
         PSV=PSV5N(S,V)
      ELSE IF (IREG .EQ. 9) THEN
         PSV=PSV9N(S,V)
      ELSE
         PSV=0.D0
      ENDIF
C
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION PSV1N(S,V)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision NUSV1N
      EXTERNAL NUSV1N
C
      EPS=1.D-9
C
      TZ1=TPS1N(60.D0,S)
      VZ1=VPT1N(60.D0,TZ1)
      TZ2=TPS1N(25.D0,S)
      VZ2=VPT1N(25.D0,TZ2)
      IF (V .LE. VZ1) THEN
         P1=60.D0
         P2=100.D0
      ELSE IF (V .LE. VZ2) THEN
         P1=25.D0
         P2=60.D0
      ELSE
         P1=PSATTN(TZ2)
         P2=25.D0
      ENDIF
C
      CALL WNSSAT3(P1,P2,NUSV1N,V,S,EPS,X,IX)
C
      PSV1N=X
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUSV1N(P,S,V)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      T=TPS1N(P,S)
      NUSV1N=VPT1N(P,T)-V
C
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION PSATS1N(S)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision NUS1SATN
      EXTERNAL NUS1SATN
C
      EPS=1.D-9
C
      IF (S .LE. 3.778281339544279D0) THEN
         P1=0.1d0
         P2=10.d0
      ELSE
         P1=18.D0
         P2=21.D0
      ENDIF
      DUM=1.D0
C
      CALL WNSSAT3(P1,P2,NUS1SATN,DUM,S,EPS,X,IX)
C
      PSATS1N=X
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUS1SATN(P,S,DUM)
C***********************************************************************
      IMPLICIT double precision (A-H,O-Z)
      COMMON/CCRIT/VC,SC,HC
C
      Z=DUM
      IF (P .GT. 22.064D0) THEN
         PN=22.064D0
      ELSE   
         PN=P
      ENDIF
      T=TSATPN(PN)
      IF (T .LE. 623.15D0) THEN
         SN=SPT1N(PN,T)
      ELSE IF (T .GE. 647.0959D0) THEN
           SN=SC
        ELSE
         CALL FSATP(DV,DL,TOUT,T,PN)
         VN=1.D0/DL
         SN=SVT3N(VN,T)
      ENDIF
      NUS1SATN=SN-S
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION SBV(V)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision SNUVBN
      EXTERNAL SNUVBN
C
      EPS=1.D-6
C
      P1=PSATTN(623.15D0)
      P2=20.D0
C
      DUM=1.d0
C
      CALL WNSSAT3(P1,P2,SNUVBN,DUM,V,EPS,X,IX)
C      
      TB=FB23P(x)
      SBV=SPT2N(X,TB)
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION SNUVBN(P,V,DUM)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      Z=DUM
      PG=PSATTN(623.15D0)
      IF (P .LT. PG) THEN
         PN=PG
      ELSE   
         PN=P
      ENDIF
      T=FB23P(PN)
      SNUVBN=VPT2N(PN,T)-V
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION PSV2N(S,V)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision NUSV2N
      EXTERNAL NUSV2N
C
      EPS=1.D-6
C
      IF (V .LT. 8.8D-3) THEN
         P1=15.D0
         P2=25.D0
      ELSEIF (V .LT. 1.8D-2) THEN
         P1=10.D0
         P2=15.D0
      ELSEIF (V .LT. 1.95D-1) THEN
         P1=1.D0
         P2=10.D0
      ELSEIF (V .LT. 1.7D0) THEN
         P1=1.D-1
         P2=1.D0
      ELSEIF (V .LT. 1.467D1) THEN
         P1=1.D-2
         P2=1.D-1
      ELSEIF (V .LT. 1.292D2) THEN
         P1=1.D-3
         P2=1.D-2
      ELSE
         P1=5.D-4
         P2=1.D-3
      ENDIF
C
      CALL WNVS3(P1,P2,NUSV2N,V,S,EPS,X,IX)
C
      PSV2N=X
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUSV2N(P,S,V)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      T=TPS2N(P,S)
      NUSV2N=VPT2N(P,T)-V
C
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION PSV3N(S,V)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision NUSV3N
      EXTERNAL NUSV3N
C
      EPS=1.D-9
C
      T1=650.d0
      T2=800.D0
C
      CALL WNVS3(T1,T2,NUSV3N,S,V,EPS,X,IX)
C
      PSV3N=PVT3N(V,X)
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUSV3N(T,V,S)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      NUSV3N=SVT3N(V,T)-S
C
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION PSV5N(S,V)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision NUSV5N
      EXTERNAL NUSV5N
C
      EPS=1.D-6
C
      IF (V .LE. 5.D0) THEN
           P1=0.1D0
           P2=10.D0
        ELSE
           P1=0.001d0
           P2=0.21d0
        ENDIF
c      P1=0.001d0
c      P2=8.D0
C
      CALL WNVS3(P1,P2,NUSV5N,V,S,EPS,X,IX)
C
      PSV5N=X
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUSV5N(P,S,V)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      T=TPS(P,S,5)
      NUSV5N=VPT5N(P,T)-V
C
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION PSV9N(S,V)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision NUSV9N
      EXTERNAL NUSV9N
      COMMON/CCRIT/VC,SC,HC
C
      EPS=1.D-6
C
      P1=0.001D0
      IF ((DABS(S-SC)) .LT. 1.d-2) THEN 
           P2=22.0d0
        ELSEIF (S .LE. SC) THEN
         P2=PSATS1N(S)
      ELSE
         P2=PSATS2N(S)
      ENDIF
C
      CALL WNVS3(P1,P2,NUSV9N,V,S,EPS,X,IX)
C
      PSV9N=X
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUSV9N(P,S,V)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      IF (P .GT. 22.064D0) THEN
         PN=22.06D0
      ELSEIF (P .LT. 5.D-4) THEN
         PN=5.D-4
      ELSE
         PN=P
      ENDIF
      T=TSATPN(PN)
      IF (T .LE. 623.15D0) THEN
         S1=SPT1N(PN,T)
         S2=SPT2N(PN,T)
         V1=VPT1N(PN,T)
         V2=VPT2N(PN,T)
      ELSE
         CALL FSATP(DV,DL,TOUT,T,PN)
         V1=1.D0/DL
         V2=1.D0/DV
         S1=SVT3N(V1,T)
         S2=SVT3N(V2,T)
      ENDIF
      NUSV9N=S1+(V-V1)/(V2-V1)*(S2-S1)-S
C
      END
C
C*******************************************************************
      SUBROUTINE WNVS3(XA,XB,F,P,T,EPS,X,IX)
C*******************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      EXTERNAL F
      X1=XA
      F1=F(X1,T,P)
      X3=XB
      F3=F(X3,T,P)
      IX=0
C++ SCHLEIFE ZUR ERMITTLUNG DER NULLSTELLE +++++++++++++++++++++++++++
      DO 100 I=1,40
      IF( F1 .NE. F3 ) THEN
       X=X1+(X3-X1)*F1/(F1-F3)
      ELSE
       IX=3
       GOTO 999
      END IF
      IF(X .LT. 0.D0) X=(X1+X3)/2.D0
      IF(DABS(X) .LT. 1.D-8) THEN
        IF(DABS(X-X1).LT.EPS) RETURN
      ELSE
        IF(DABS((X-X1)/X).LT.EPS) RETURN
      END IF
      F2=F(X,T,P)
      X2=X1-(X1-X3)/2.D0
      P1=F2*F1
      P3=F2*F3
      IF ((P1 .LT. 0.D0) .AND. (P3 .LT. 0.D0)) THEN
        P1=DABS(P1)
        P3=DABS(P3)
      ENDIF
      IF(P1 .LE. P3) THEN
        X3=X1
        F3=F1
      END IF
      X1=X
      F1=F2
      IF((X2-X3)*(X2-X1) .GE. 0.D0) GOTO 100
c      X=(X1+X3)/2.D0
      IF ((DABS(F1/F3)) .LE. (0.3D0)) THEN
         X=X1+DABS(F1/F3)*1.5D0*(X3-X1)
         IREM=1
      ELSE IF ((DABS(F1/F3)) .GE. (3.D0)) THEN
         IREM=2
      ELSE
         X=(X1+X3)/2.D0
         IREM=0
      ENDIF
      F2=F(X,T,P)
      IF ((((F2*F1) .GE. 0.D0) .AND. (IREM .EQ. 1)) .OR.
     1    (((F2*F3) .GE. 0.D0) .AND. (IREM .EQ. 2))) THEN
         X=(X1+X3)/2.D0
         F2=F(X,T,P)
      ENDIF
      IF(F2*F1 .LE. F2*F3) THEN
        X3=X1
        F3=F1
      END IF
      X1=X
      F1=F2
 100  CONTINUE
C++ ENDE DER SCHLEIFE ++++++++++++++++++++++++++++++++++++++++++++++++
      IX=1
 999  RETURN
      END
C
C*******************************************************************
      SUBROUTINE WNSSAT3(XA,XB,F,P,T,EPS,X,IX)
C*******************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      EXTERNAL F
      X1=XA
      F1=F(X1,T,P)
      X3=XB
      F3=F(X3,T,P)
      IX=0
C++ SCHLEIFE ZUR ERMITTLUNG DER NULLSTELLE +++++++++++++++++++++++++++
      DO 100 I=1,40
      IF( F1 .NE. F3 ) THEN
       X=X1+(X3-X1)*F1/(F1-F3)
      ELSE
       IX=3
       GOTO 999
      END IF
      IF(X .LT. 0.D0) X=(X1+X3)/2.D0
      IF(DABS(X) .LT. 1.D-8) THEN
        IF(DABS(X-X1).LT.EPS) RETURN
      ELSE
        IF(DABS((X-X1)/X).LT.EPS) RETURN
      END IF
      F2=F(X,T,P)
      X2=X1-(X1-X3)/2.D0
      P1=F2*F1
      P3=F2*F3
      IF ((P1 .LT. 0.D0) .AND. (P3 .LT. 0.D0)) THEN
        P1=DABS(P1)
        P3=DABS(P3)
      ENDIF
      IF(P1 .LE. P3) THEN
        X3=X1
        F3=F1
      END IF
      X1=X
      F1=F2
      IF((X2-X3)*(X2-X1) .GE. 0.D0) GOTO 100
c      X=(X1+X3)/2.D0
      IF ((DABS(F1/F3)) .LE. (0.3D0)) THEN
         X=X1+DABS(F1/F3)*1.5D0*(X3-X1)
         IREM=1
      ELSE IF ((DABS(F1/F3)) .GE. (3.D0)) THEN
         IREM=2
      ELSE
         X=(X1+X3)/2.D0
         IREM=0
      ENDIF
      F2=F(X,T,P)
      IF ((((F2*F1) .GE. 0.D0) .AND. (IREM .EQ. 1)) .OR.
     1    (((F2*F3) .GE. 0.D0) .AND. (IREM .EQ. 2))) THEN
         X=(X1+X3)/2.D0
         F2=F(X,T,P)
      ENDIF
      IF(F2*F1 .LE. F2*F3) THEN
        X3=X1
        F3=F1
      END IF
      X1=X
      F1=F2
 100  CONTINUE
C++ ENDE DER SCHLEIFE ++++++++++++++++++++++++++++++++++++++++++++++++
      IX=1
 999  RETURN
      END
C
C***********************************************************************
      SUBROUTINE REGSHS(H,S,IREG)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      COMMON/CCRIT/VC,SC,HC
C
      double precision HOLD, SOLD
      SAVE HOLD, SOLD
      INTEGER IREGOLD
      SAVE IREGOLD
      data HOLD, SOLD, IREGOLD / -1.D-6, -1.D-6, -1 /
C
      IF ((DABS(H-HOLD).LT. 1.D-6) .AND.
     *   (DABS(S-SOLD).LT. 1.D-6)) THEN
         IREG=IREGOLD
         GOTO 999
      END IF
C
      IF ((H .LT. -0.05D0) .OR. (S .LT. -0.5D0) .OR.
     &    (H .GE. 7377.D0) .OR. (S .GT. 13.99764756652987D0)) THEN
         IREG=0
         GOTO 999
      ENDIF
C
        IREG=0
C
      IF ((H .LT. HC) .AND. (S .LT. SC)) THEN
         PS=PSATH1N(H)
         TS=TSATPN(PS)
         IF (TS .LE. 623.15D0) THEN
            S1=SPT1N(PS,TS)
         ELSE
            CALL FSATP(DV,DL,TOUT,TS,PS)
            V1=1.D0/DL
            S1=SVT3N(V1,TS)
         ENDIF   
         IF (S .GT. S1) THEN
            IREG=9
         ELSEIF (H .LE. 1553.922503371486D0) THEN
            TN=TPH1N(100.D0,H)
            SN=SPT1N(102.D0,TN)
            IF (S .GE. SN) THEN
               IREG=1
            ELSE
               IREG=0
            ENDIF
         ELSEIF (H .GE. 1670.858218274385D0) THEN
            TN=TPH3N(100.D0,H)
            VN=VPT3N(102.D0,TN)
            SN=SVT3N(VN,TN)
            IF (S .GE. SN) THEN
               IREG=3
            ELSE
               IREG=0
            ENDIF
         ELSE
ccc        Grenze 1/3
C      
         HG1 = (((((1.2985145661014D-3 * S - 2.22001959072D-2) * S
     &     +1.14738086963407D-1) * S - 364.83398463958752D0) * S
     &     +4568.3374037778948195D0) * S - 18377.982160634557656D0) * S
     &     +25560.70255570028708D0
C
            IF (H .LT. HG1) THEN
               TN=TPH1N(100.D0,H)
               SN=SPT1N(102.D0,TN)
               IF (S .GE. SN) THEN
                  IREG=1
               ELSE
                  IREG=0
               ENDIF
            ELSE
               HG3=HG1+3.d0
               PHF=PHS3N(H,S)
               IF (H .GT. HG3) THEN
                            IF (PHF .LE. 102.D0) THEN
                     IREG=3
                  ELSE
                     IREG=0
                  ENDIF
               ELSE
                  THF=TPS3N(PHF,S)
                  IF ((THF .LE. 623.15D0) .AND. (PHF .LE. 102.D0)) THEN
                     IREG=1
                  ELSEIF (PHF .LE. 102.D0) THEN
                     IREG=3
                  ELSE
                     IREG=0
                    ENDIF
               ENDIF
            ENDIF  
         ENDIF
ccc
      ELSEIF ((H .GE. HC) .AND. (S .LT. 5.0465D0)) THEN
         IF (S .LT. 4.185895623423D0) THEN
            IREG=0
         ELSEIF (S .LE. SC) THEN         
            TN=TPH3N(100.D0,H)
            VN=VPT3N(102.D0,TN)
            SN=SVT3N(VN,TN)
            IF (S .GE. SN) THEN
               IREG=3
            ELSE   
               IREG=0
            ENDIF
         ELSE
            PS=PSATS2N(S)
            TS=TSATPN(PS)
            IF (TS .LE. 623.15D0) THEN
               H2=HPT2N(PS,TS)
            ELSE
               CALL FSATP(DV,DL,TOUT,TS,PS)
               V2=1.D0/DV
               H2=HVT3N(V2,TS)
            ENDIF
            IF (H .LE. H2) THEN
               IREG=9
            ELSE
               TN=TPH3N(100.D0,H)
               VN=VPT3N(102.D0,TN)
               SN=SVT3N(VN,TN)
               IF (S .GE. SN) THEN
                  IREG=3
               ELSE   
                  IREG=0
               ENDIF
            ENDIF
         ENDIF
      ELSE
         IF (S .LE. 9.15549147374441d0) THEN
            PS=PSATS2N(S)
            TS=TSATPN(PS)
              IF (TS .LE. 623.15D0) THEN
               H2=HPT2N(PS,TS)
            ELSE
               CALL FSATP(DV,DL,TOUT,TS,PS)
               V2=1.D0/DV
               H2=HVT3N(V2,TS)
            ENDIF
           ELSE
              H2=HC
           ENDIF
C
         IF (H .LE. H2) THEN
            IREG=9
         ELSEIF ((S .LE. 5.2615D0) .AND. (H .LE. 2820.D0)) THEN
ccc      Grenze 2/3 
            PN2=PHS2N(H,S)
            TB=FB23P(PN2)
            TN2=TPS2N(PN2,S)
            IF (TN2 .GT. TB) THEN
               IF (PN2 .LE. 102.D0) THEN
                  IREG=2
               ELSE
                  IREG=0
               ENDIF
            ELSE
               PN3=PHS3N(H,S)
               IF (PN3 .LE. 102.D0) THEN
                  IREG=3
               ELSE
                  IREG=0
               ENDIF
            ENDIF
         ELSEIF ((S .LE. 6.040483671712074D0) .AND. 
     &    (H .LT. 3750.D0)) THEN
            TN=TPH2N(100.D0,H)
            SN=SPT2N(102.D0,TN)
            IF (S .GE. SN) THEN
               IREG=2
            ELSE   
               IREG=0
            ENDIF
         ELSE
            H100=3715.18894d0
            H90=3753.01631d0
            H80=3793.32250d0
            H70=3835.81423d0
            H60=3880.15394d0
            H50=3925.96041d0
            H40=3972.80939d0
            H30=4020.23405d0
            H20=4067.72544d0
            H10=4114.73278d0
            HD0=4156.13678d0
            HD1=4160.21176d0
            HD2=4160.61850d0
            HD3=4160.65917d0
            HD4=4160.66324d0
            S100=6.04048367d0
            S90=6.11837329d0
            S80=6.20385037d0
            S70=6.29820442d0
            S60=6.40340896d0
            S50=6.52264231d0
            S40=6.66140160d0
            S30=6.83025450d0
            S20=7.05335918d0
            S10=7.40867489d0
            SD0=8.50236101d0
            SD1=9.56810070d0
            SD2=10.6311066d0
            SD3=11.6938398d0
            SD4=12.7565457d0
            IF (H .LT. H100) THEN
               IREG=2
            ELSEIF ((H .GE. H100) .AND. (H .LE. H90)) THEN
               IF (S .GE. S90) THEN
                  IREG=2
               ELSEIF (S .LE. S100) THEN
                  IREG=0
               ELSE
                  PN=PHS2N(H,S)
                  TN=TPS2N(PN,S)
                  IF (TN .LE. 1073.15D0) THEN
                     IREG=2
                  ELSE
                     IREG=0
                  ENDIF
               ENDIF
            ELSEIF ((H .GT. H90) .AND. (H .LE. H80)) THEN
               IF (S .GE. S80) THEN
                  IREG=2
               ELSEIF (S .LE. S90) THEN
                  IREG=0
               ELSE
                  PN=PHS2N(H,S)
                  TN=TPS2N(PN,S)
                  IF (TN .LE. 1073.15D0) THEN
                     IREG=2
                  ELSE
                     IREG=0
                  ENDIF
               ENDIF
            ELSEIF ((H .GT. H80) .AND. (H .LE. H70)) THEN
               IF (S .GE. S70) THEN
                  IREG=2
               ELSEIF (S .LE. S80) THEN
                  IREG=0
               ELSE
                  PN=PHS2N(H,S)
                  TN=TPS2N(PN,S)
                  IF (TN .LE. 1073.15D0) THEN
                     IREG=2
                  ELSE
                     IREG=0
                  ENDIF
               ENDIF
            ELSEIF ((H .GT. H70) .AND. (H .LE. H60)) THEN
               IF (S .GE. S60) THEN
                  IREG=2
               ELSEIF (S .LE. S70) THEN
                  IREG=0
               ELSE
                  PN=PHS2N(H,S)
                  TN=TPS2N(PN,S)
                  IF (TN .LE. 1073.15D0) THEN
                     IREG=2
                  ELSE
                     IREG=0
                  ENDIF
               ENDIF
            ELSEIF ((H .GT. H60) .AND. (H .LE. H50)) THEN
               IF (S .GE. S50) THEN
                  IREG=2
               ELSEIF (S .LE. S60) THEN
                  IREG=0
               ELSE
                  PN=PHS2N(H,S)
                  TN=TPS2N(PN,S)
                  IF (TN .LE. 1073.15D0) THEN
                     IREG=2
                  ELSE
                     IREG=0
                  ENDIF
               ENDIF
            ELSEIF ((H .GT. H50) .AND. (H .LE. H40)) THEN
               IF (S .GE. S40) THEN
                  IREG=2
               ELSEIF (S .LE. S50) THEN
                  IREG=0
               ELSE
                  PN=PHS2N(H,S)
                  TN=TPS2N(PN,S)
                  IF (TN .LE. 1073.15D0) THEN
                     IREG=2
                  ELSE
                     IREG=0
                  ENDIF
               ENDIF
            ELSEIF ((H .GT. H40) .AND. (H .LE. H30)) THEN
               IF (S .GE. S30) THEN
                  IREG=2
               ELSEIF (S .LE. S40) THEN
                  IREG=0
               ELSE
                  PN=PHS2N(H,S)
                  TN=TPS2N(PN,S)
                  IF (TN .LE. 1073.15D0) THEN
                     IREG=2
                  ELSE
                     IREG=0
                  ENDIF
               ENDIF
            ELSEIF ((H .GT. H30) .AND. (H .LE. H20)) THEN
               IF (S .GE. S20) THEN
                  IREG=2
               ELSEIF (S .LE. S30) THEN
                  IREG=0
               ELSE
                  PN=PHS2N(H,S)
                  TN=TPS2N(PN,S)
                  IF (TN .LE. 1073.15D0) THEN
                     IREG=2
                  ELSE
                     IREG=0
                  ENDIF
               ENDIF
            ELSEIF ((H .GT. H20) .AND. (H .LE. H10)) THEN
               IF (S .GE. S10) THEN
                  IREG=2
               ELSEIF (S .LE. S20) THEN
                  IREG=0
               ELSE
                  PN=PHS2N(H,S)
                  TN=TPS2N(PN,S)
                  IF (TN .LE. 1073.15D0) THEN
                     IREG=2
                  ELSE
                     IREG=0
                  ENDIF
               ENDIF
            ELSEIF ((H .GT. H10) .AND. (H .LE. HD0)) THEN
               IF (S .GE. SD0) THEN
                  IREG=2
               ELSEIF (S .LE. S10) THEN
ccccc            
                  TN=TPH5N(10.D0,H)
                  SN=SPT5N(10.2D0,TN)
                  IF (S .GE. SN) THEN
                     IREG=5
                  ELSE
                     IREG=0
                  ENDIF
cccc
               ELSE
                  PN=PHS2N(H,S)
                  TN=TPS2N(PN,S)
                  IF (TN .LE. 1073.15D0) THEN
                     IREG=2
                  ELSE
                     TN=TPH5N(10.D0,H)
                     SN=SPT5N(10.2D0,TN)
                     IF (S .GE. SN) THEN
                        IREG=5
                     ELSE
                        IREG=0
                     ENDIF
                  ENDIF
               ENDIF
            ELSEIF ((H .GT. HD0) .AND. (H .LE. HD1)) THEN
               IF (S .GE. SD1) THEN
                  IREG=2
               ELSEIF (S .LE. SD0) THEN
ccccc            
                  TN=TPH5N(10.D0,H)
                  SN=SPT5N(10.2D0,TN)
                  IF (S .GE. SN) THEN
                     IREG=5
                  ELSE
                     IREG=0
                  ENDIF
cccc
               ELSE
                  PN=PHS2N(H,S)
                  TN=TPS2N(PN,S)
                  IF (TN .LE. 1073.15D0) THEN
                     IREG=2
                  ELSE
                     TN=TPH5N(10.D0,H)
                     SN=SPT5N(10.2D0,TN)
                     IF (S .GE. SN) THEN
                        IREG=5
                     ELSE
                        IREG=0
                     ENDIF
                  ENDIF
               ENDIF
            ELSEIF ((H .GT. HD1) .AND. (H .LE. HD2)) THEN
               IF (S .GE. SD2) THEN
                  IREG=2
               ELSEIF (S .LE. SD1) THEN
ccccc            
                  TN=TPH5N(10.D0,H)
                  SN=SPT5N(10.2D0,TN)
                  IF (S .GE. SN) THEN
                     IREG=5
                  ELSE
                     IREG=0
                  ENDIF
cccc
               ELSE
                  PN=PHS2N(H,S)
                  TN=TPS2N(PN,S)
                  IF (TN .LE. 1073.15D0) THEN
                     IREG=2
                  ELSE
                     TN=TPH5N(10.D0,H)
                     SN=SPT5N(10.2D0,TN)
                     IF (S .GE. SN) THEN
                        IREG=5
                     ELSE
                        IREG=0
                     ENDIF
                  ENDIF
               ENDIF
            ELSEIF ((H .GT. HD2) .AND. (H .LE. HD3)) THEN
               IF (S .GE. SD3) THEN
                  IREG=2
               ELSEIF (S .LE. SD2) THEN
ccccc            
                  TN=TPH5N(10.D0,H)
                  SN=SPT5N(10.2D0,TN)
                  IF (S .GE. SN) THEN
                     IREG=5
                  ELSE
                     IREG=0
                  ENDIF
cccc
               ELSE
                  PN=PHS2N(H,S)
                  TN=TPS2N(PN,S)
                  IF (TN .LE. 1073.15D0) THEN
                     IREG=2
                  ELSE
                     TN=TPH5N(10.D0,H)
                     SN=SPT5N(10.2D0,TN)
                     IF (S .GE. SN) THEN
                        IREG=5
                     ELSE
                        IREG=0
                     ENDIF
                  ENDIF
               ENDIF
            ELSEIF ((H .GT. HD3) .AND. (H .LE. HD4)) THEN
               IF (S .GE. SD4) THEN
                  IREG=2
               ELSEIF (S .LE. SD3) THEN
ccccc            
                  TN=TPH5N(10.D0,H)
                  SN=SPT5N(10.2D0,TN)
                  IF (S .GE. SN) THEN
                     IREG=5
                  ELSE
                     IREG=0
                  ENDIF
cccc
               ELSE
                  PN=PHS2N(H,S)
                  TN=TPS2N(PN,S)
                  IF (TN .LE. 1073.15D0) THEN
                     IREG=2
                  ELSE
                     TN=TPH5N(10.D0,H)
                     SN=SPT5N(10.2D0,TN)
                     IF (S .GE. SN) THEN
                        IREG=5
                     ELSE
                        IREG=0
                     ENDIF
                  ENDIF
               ENDIF
            ELSEIF (H .LE. 7374.751711911952D0) THEN
               TN=TPH5N(10.D0,H)
               SN=SPT5N(10.2D0,TN)
               IF (S .GE. SN) THEN
                  IREG=5
               ELSE
                  IREG=0
               ENDIF
ccc            
            ELSEIF (H .LE. 7376.540702013396D0) THEN
               PN=PHS5N(H,S)
               TN=TPH5N(PN,H)
               IF (TN .LE. 2273.2D0) THEN
                  IREG=5
               ELSE
                  IREG=0
               ENDIF
            ELSE
               IREG=0
            ENDIF
         ENDIF
      ENDIF
C
      IF (IREG .EQ. 9) THEN
         TS=TSATPN(5.D-4)
         S1=SPT1N(5.D-4,TS)
         S2=SPT2N(5.D-4,TS)
         H1=HPT1N(5.D-4,TS)
         H2=HPT2N(5.D-4,TS)
         XS=(S-S1)/(S2-S1)
         XH=(H-H1)/(H2-H1)
         IF (XH .LT. XS) THEN
            IREG=0
         ENDIF
      ENDIF
C
      IF ((IREG .EQ. 2) .AND. (S .GT. 9.1D0)) THEN
         PV=PHS2N(H,S)
         IF (PV .LT. 5.D-4) THEN
            IREG=0
         ENDIF
      ENDIF
C
999   CONTINUE   
C
      HOLD=H
      SOLD=S
      IREGOLD=IREG
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION PHS(H,S,IREG)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      IF (IREG .EQ. 0) THEN
         CALL REGSHS(H,S,IREG)
      ENDIF
C  
      IF (IREG .EQ. 1) THEN
         PHS=PHS1N(H,S)
      ELSE IF (IREG .EQ. 2) THEN
         PHS=PHS2N(H,S)
      ELSE IF (IREG .EQ. 3) THEN
         PHS=PHS3N(H,S)
      ELSE IF (IREG .EQ. 9) THEN
         PHS=PHS9N(H,S)
      ELSE IF (IREG .EQ. 5) THEN
         PHS=PHS5N(H,S)
      ELSE
         PHS=0.D0
      ENDIF
C
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION PHS1N(H,S)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      PARAMETER(a01=-.691997014660582D+00)
      PARAMETER(a02=-.183612548787560D+02)
      PARAMETER(a03=-.928332409297335D+01)
      PARAMETER(a04= .659639569909906D+02)
      PARAMETER(a05=-.162060388912024D+02)
      PARAMETER(a06= .450620017338667D+03)
      PARAMETER(a07= .854680678224170D+03)
      PARAMETER(a08= .607523214001162D+04)
      PARAMETER(a09= .326487682621856D+02)
      PARAMETER(a10=-.269408844582931D+02)
      PARAMETER(a11=-.319947848334300D+03)
      PARAMETER(a12=-.928354307043320D+03)
      PARAMETER(a13= .303634537455249D+02)
      PARAMETER(a14=-.650540422444146D+02)
      PARAMETER(a15=-.430991316516130D+04)
      PARAMETER(a16=-.747512324096068D+03)
      PARAMETER(a17= .730000345529245D+03)
      PARAMETER(a18= .114284032569021D+04)
      PARAMETER(a19=-.436407041874559D+03)
C
      h1= h*(1.d0/3400.d0)+0.05d0
      h2= h1*h1
      h3= h1*h2
      s1= s*(1.d0/7.6d0)+0.05d0
      s2= s1*s1
      s4= s2*s2
c
       PHS1n= (               a01
     *      +        h1*     (a09
     *      +        h1*     (a13
     *      +        h3*      a19))
     *      +s1*(            (a02
     *      +        h1*     (a10
     *      +        h1*     (a14
     *      +        h2*      a17)))
     *      +s1*(             a03
     *      +s2*(             a04
     *      +        h1*     (a11
     *      +        h2*     (a16
     *      +        h1*      a18))
     *      +s1*(             a05
     *      +s1*(            (a06
     *      +        h1*      a12)
     *      +s2*(             a07
     *      +s2*(    h2*      a15
     *      +s4*              a08))))))))*100.d0
C
        END
C************************************************************************
      DOUBLE PRECISION FUNCTION PHS2N(H,S)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      PARAMETER (a01=-.182575361923032D-01)
      PARAMETER (a02=-.125229548799536D+00)    
      PARAMETER (a03= .592290437320145D+00)    
      PARAMETER (a04= .604769706185122D+01)    
      PARAMETER (a05= .238624965444474D+03)    
      PARAMETER (a06=-.298639090222922D+03)    
      PARAMETER (a07= .512250813040750D-01)
      PARAMETER (a08=-.437266515606486D+00)    
      PARAMETER (a09= .413336902999504D+00)    
      PARAMETER (a10=-.516468254574773D+01)    
      PARAMETER (a11=-.557014838445711D+01)    
      PARAMETER (a12= .128555037824478D+02)    
      PARAMETER (a13= .114144108953290D+02)    
      PARAMETER (a14=-.119504225652714D+03)    
      PARAMETER (a15=-.284777985961560D+04)    
      PARAMETER (a16= .431757846408006D+04)    
      PARAMETER (a17= .112894040802650D+01)    
      PARAMETER (a18= .197409186206319D+04)    
      PARAMETER (a19= .151612444706087D+04)    
      PARAMETER (a20= .141324451421235D-01)
      PARAMETER (a21= .585501282219601D+00)    
      PARAMETER (a22=-.297258075863012D+01)    
      PARAMETER (a23= .594567314847319D+01)    
      PARAMETER (a24=-.623656565798905D+04)    
      PARAMETER (a25= .965986235133332D+04)    
      PARAMETER (a26= .681500934948134D+01)    
      PARAMETER (a27=-.633207286824489D+04)    
      PARAMETER (a28=-.558919224465760D+01)    
      PARAMETER (a29= .400645798472063D-01)

      PARAMETER (b01= .801496989929495D-01)
      PARAMETER (b02=-.543862807146111D+00)    
      PARAMETER (b03= .337455597421283D+00)    
      PARAMETER (b04= .890555451157450D+01)    
      PARAMETER (b05= .313840736431485D+03)    
      PARAMETER (b06= .797367065977789D+00)    
      PARAMETER (b07=-.121616973556240D+01)    
      PARAMETER (b08= .872803386937477D+01)    
      PARAMETER (b09=-.169769781757602D+02)    
      PARAMETER (b10=-.186552827328416D+03)    
      PARAMETER (b11= .951159274344237D+05)    
      PARAMETER (b12=-.189168510120494D+02)    
      PARAMETER (b13=-.433407037194840D+04)    
      PARAMETER (b14= .543212633012715D+09)    
      PARAMETER (b15= .144793408386013D+00)    
      PARAMETER (b16= .128024559637516D+03)    
      PARAMETER (b17=-.672309534071268D+05)    
      PARAMETER (b18= .336972380095287D+08)    
      PARAMETER (b19=-.586634196762720D+03)    
      PARAMETER (b20=-.221403224769889D+11)    
      PARAMETER (b21= .171606668708389D+04)    
      PARAMETER (b22=-.570817595806302D+09)    
      PARAMETER (b23=-.312109693178482D+04)    
      PARAMETER (b24=-.207841384633010D+07)    
      PARAMETER (b25= .305605946157786D+13)    
      PARAMETER (b26= .322157004314333D+04)    
      PARAMETER (b27= .326810259797295D+12)    
      PARAMETER (b28=-.144104158934487D+04)    
      PARAMETER (b29= .410694867802691D+03)    
      PARAMETER (b30= .109077066873024D+12)    
      PARAMETER (b31=-.247964654258893D+14)    
      PARAMETER (b32= .188801906865134D+10)    
      PARAMETER (b33=-.123651009018773D+15)    
c
      PARAMETER (c01= .112225607199012D+00)     
      PARAMETER (c02=-.339005953606712d+01)     
      PARAMETER (c03=-.320503911730094d+02)     
      PARAMETER (c04=-.197597305104900d+03)     
      PARAMETER (c05=-.407693861553446d+03)     
      PARAMETER (c06= .132943775222331d+05)     
      PARAMETER (c07= .170846839774007d+01)     
      PARAMETER (c08= .373694198142245d+02)     
      PARAMETER (c09= .358144365815434d+04)     
      PARAMETER (c10= .423014446424664d+06)     
      PARAMETER (c11=-.751071025760063d+09)     
      PARAMETER (c12= .523446127607898d+02)     
      PARAMETER (c13=-.228351290812417d+03)     
      PARAMETER (c14=-.960652417056937d+06)     
      PARAMETER (c15=-.807059292526074d+08)     
      PARAMETER (c16= .162698017225669d+13)     
      PARAMETER (c17= .772465073604171d+00)     
      PARAMETER (c18= .463929973837746d+05)     
      PARAMETER (c19=-.137317885134128d+08)     
      PARAMETER (c20= .170470392630512d+13)     
      PARAMETER (c21=-.251104628187308d+14)     
      PARAMETER (c22= .317748830835520d+14)     
      PARAMETER (c23= .538685623675312d+02)     
      PARAMETER (c24=-.553089094625169d+05)     
      PARAMETER (c25=-.102861522421405d+07)     
      PARAMETER (c26= .204249418756234d+13)     
      PARAMETER (c27= .273918446626977d+09)     
      PARAMETER (c28=-.263963146312685d+16) 
      PARAMETER (c29=-.107890854108088d+10)     
      PARAMETER (c30=-.296492620980124d+11)     
      PARAMETER (c31=-.111754907323424d+16) 
C
      PARAMETER (d1=-3498.98083432139d0)
      PARAMETER (d2=2575.60716905876d0)
      PARAMETER (d3=-421.073558227969d0)
      PARAMETER (d4=27.6349063799944d0)
c
      PARAMETER (hq4200= 1.d0/4200.d0,
     1           sq12= 1.d0/12.d0)
C
      PARAMETER (hq4100= 1.D0/4100.D0,
     1           sq79= 1.d0/7.9d0)
C
      PARAMETER (hq3500= 1.d0/3500.d0,
     1           sq59= 1.d0/5.9d0)       
C
        IF (s.LT.5.85d0) GOTO 30
        IF (s.LT.6.06970915951902d0) GOTO 20
        IF (s.GE.7.85234039987821d0) GOTO 10      
        h4mp= d1+s*(d2+s*(d3+s*d4))
        IF (h.GT.h4mp) GOTO 20
C
C Subregion a
C
10      h1= h*hq4200-.5d0
        h2= h1*h1
      h6= h2*h2*h2
        s1= s*sq12-1.2d0
        s2= s1*s1
        s4= s2*s2
        s6= s2*s4
        phs=         h1*     (a07
     *      +        h2*      a20)
     *      +s1*             (a01
     *      +        h1*     (a08
     *      +        h6*      a29))
     *      +s2*(    h1*     (a09
     *      +        h2*      a21)
     *      +s1*(            (a02
     *      +        h1*     (a10
     *      +        h1*     (a17
     *      +        h1*     (a22         
     *      +        h2*     (a26
     *      +        h1*      a28)))))
     *      +s2*     h1*      a11))
     *      +s6*(            (a03
     *      +        h1*     (a12
     *      +        h2*      a23))
     *      +s4*(    h1*      a13
     *      +s6*(            (a04
     *      +        h1*     (a14
     *      +        h1*     (a18
     *      +        h1*     (a24
     *      +        h1*     (a25
     *      +        h1*      a27)))))
     *      +s4*(            (a05
     *      +        h1*     (a15
     *      +        h1*      a19))
     *      +s2*             (a06
     *      +        h1*                a16)))))
C
        phs2= phs*phs
      PHS2n= phs2*phs2*4.d0
        GOTO 99
C
C Subregion 2b
C
20      h1= h*hq4100-.6d0
        h2= h1*h1
        h3= h1*h2
        h4= h1*h3
        h6= h3*h3
        h7= h1*h6
        h12= h6*h6
        s1= s*sq79-1.01d0
        s2=  s1*s1
        phs=                  b01
     *      +        h1*     (b06
     *      +        h2*      b15)
     *      +s1*   (          b02
     *      +        h1*     (b07
     *      +        h1*     (b12
     *      +        h1*     (b16
     *      +        h1*     (b19 
     *      +        h1*     (b21 
     *      +        h1*     (b23 
     *      +        h1*     (b26 
     *      +        h1*      b28))))))) 
     *      +s1*   (          b03
     *      +        h1*      b08
     *      +s1*   ( h1*     (b09
     *      +        h7*      b29)
     *      +s1*   (          b04
     *      +s1*   ( h1*      b10
     *      +s1*   ( h2*      b13
     *      +s1*   ( h3*      b17
     *      +s1*   (          b05
     *      +        h6*      b24
     *      +s2*   ( h12*     b32
     *      +s2*   ( h1*     (b11
     *      +        h2*     (b18
     *      +        h2*      b22)) 
     *      +s2*   ( h1*h7*   b30
     *      +s2*   ( h4*     (b20
     *      +        h3*     (b27
     *      +        h7*      b33))
     *      +s2*   ( h2*     (b14
     *      +        h4*     (b25  
     *      +              h2*      b31)))))))))))))))
C     
      phs2 = phs*phs
      PHS2n = phs2*phs2*100.d0
      GOTO 99
C
C subregion 2c
C
30      h1= h*hq3500-.7d0
        h2= h1*h1
        h3= h1*h2
        h4= h2*h2
        h5= h1*h4
        h8= h4*h4
        h14= h4*h5*h5
        s1= s*sq59-1.1d0
        s2= s1*s1
        s4= s2*s2
        phs=                  c01
     *      +        h1*     (c07
     *      +        h2*      c17)
     *      +s1*(             c02
     *      +        h5*      c23
     *      +s1*(             c03
     *      +        h1*     (c08
     *      +        h1*      c12)
     *      +s1*(             c04 
     *      +        h2*      c13
     *      +s1*(             c05
     *      +        h5*      c24
     *      +s1*(    h1*     (c09
     *      +        h2*      c18)
     *      +s1*(    h5*      c25
     *      +s1*(    h2*     (c14
     *      +        h8*     (c29
     *      +        h2*      c30))
     *      +s1*(             c06
     *      +        h1*     (c10
     *      +        h2*     (c19
     *      +        h3*      c27))
     *      +s2*(    h2*     (c15
     *      +        h14*     c31)
     *      +s4*(    h1*     (c11
     *      +        h4*      c26)
     *      +s2*(    h3*      c20
     *      +s2*(    h2*     (c16
     *      +        h1*     (c21
     *      +        h1*     (c22     
     *      +        h2*      c28)))))))))))))))
      phs2= phs*phs
      PHS2n= phs2*phs2*100.d0
99    CONTINUE
C
      END
C***********************************************************************
        DOUBLEPRECISION FUNCTION PHS3N(H,S)
C***********************************************************************
      IMPLICIT double precision (A-H,O-Z)
C
      double precision NUHS3N
      EXTERNAL NUHS3N
C
        parameter(hc=2087.54684511650d0)
        parameter(sc=4.41202148223476d0)
        parameter(a1=  651.005569610476D0  )
        parameter(a2= -2362.95383657093D0  )
        parameter(a3= -12997.2339318628D0  )
        parameter(a4=  16492.4543343504D0  )
        parameter(a5=  72192133.9249612D0  )
        parameter(a6=  1749.80414537544D0  )
        parameter(a7=  14792.0068422084D0  )
        parameter(a8= -10149.1382284543D0  )
        parameter(a9= -4103.29323474447D0  )
        parameter(b1= 0.305267197743091D-02)
        parameter(b2= 0.645158376417871D-01)
        parameter(b3=  1.12827412377111D0  )
        parameter(b4=  7.91743842409206D0  )
        parameter(b5=-0.365508795584734D-01)
        parameter(b6= -1.38502121114191D0  )
        parameter(b7= -15.0249237969224D0  )
        parameter(b8= 0.430700230780091D0  )
        parameter(b9=  9.53201719507785D0  )
        parameter(b0= -2.02182898091898D0  )
C
      EPS=1.D-9
C
        hr=h*(1.d0/hc)-1.d0
        sr=s*(1.d0/sc)-1.d0
        s2=sr*sr
C  
      TE=a1+sr*(a2+sr*(a3+sr*(a4+s2*s2*s2*sr*a5)))
     +  +hr*(a6+sr*(a7+sr*a8)+hr*a9)
C
      VE=b1+sr*(b2+sr*(b3+sr*b4))
     +  +hr*(b5+sr*(b6+sr*b7)+hr*(b8+sr*b9+hr*b0))
C
      PE=PVT3N(VE,TE)
      P1=PE*0.98D0
      P2=PE*1.02D0
C
      CALL WNHS3(P1,P2,NUHS3N,S,H,EPS,X,IX)
C      
      PHS3N=X
C  
      END
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUHS3N(P,H,S)
C***********************************************************************
      IMPLICIT double precision (A-H,O-Z)
C
      TN=TPH3N(P,H)
      VN=VPT3N(P,TN)
      NUHS3N=SVT3N(VN,TN)-S
C      
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION PHS5N(H,S)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision NUHS5N
      EXTERNAL NUHS5N
C
      EPS=1.D-9
C
      IF (S .LE. 8.D0) THEN
         P1=3.D0
         P2=10.D0
      ELSE IF (S .LE. 9.D0) THEN
         P1=0.1D0
         P2=10.D0
      ELSE IF (S .LE. 10.D0) THEN
         P1=5.D-2
         P2=10.D0
      ELSEIF (S .LE. 11.D0) THEN
         P1=5.D-3
         P2=5.D0
      ELSE IF (S .LE. 12.D0) THEN
         P1=5.D-4
         P2=1.D0
      ELSE
         P1=5.D-4
         P2=0.1D0
      ENDIF
C
      CALL WNHS3(P1,P2,NUHS5N,S,H,EPS,X,IX)
C
      IF ((X .GE. 5.D-4) .AND. (X .LE. 10.D0)) THEN
         PHS5N=X
      ELSEIF (X .LE. 10.5D0) THEN
         PHS5N=10.D0
      ELSE
         PHS5N=0.D0
      ENDIF
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUHS5N(P,H,S)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      IF (P .GT. 10.D0) THEN
         PN=10.D0
      ELSEIF (P .LT. 5.D-4) THEN
         PN=5.D-4
      ELSE
         PN=P
      ENDIF
      T=TPH5N(PN,H)
      NUHS5N=SPT5N(PN,T)-S
C
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION PHS9N(H,S)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision NUHS9N
      EXTERNAL NUHS9N
      COMMON/CCRIT/VC,SC,HC
C
      EPS=1.D-6
C
      IF (H .LE. HC) THEN
         PG=PSATH1N(H)
      ELSE
         PG=PSATS2N(S)
      ENDIF
      IF (PG .GE. 1.D-3) THEN
           P1=1.D-3
           P2=PG
        ELSE
           P1=6.11657D-4
           P2=PG
        ENDIF
C
      CALL WNHS3(P1,P2,NUHS9N,S,H,EPS,X,IX)
C
      PHS9N=X
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUHS9N(P,H,S)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      IF (P .GT. 22.064D0) THEN
         PN=22.06D0
      ELSEIF (P .LT. 5.D-4) THEN
         PN=5.D-4
      ELSE
         PN=P
      ENDIF
      T=TSATPN(PN)
      IF (T .LE. 623.15D0) THEN
         H1=HPT1N(PN,T)
         H2=HPT2N(PN,T)
         S1=SPT1N(PN,T)
         S2=SPT2N(PN,T)
        tout=t
      ELSE
         CALL FSATP(DV,DL,TOUT,T,PN)
         V1=1.D0/DL
         V2=1.D0/DV
         H1=HVT3N(V1,T)
         H2=HVT3N(V2,T)
         S1=SVT3N(V1,T)
         S2=SVT3N(V2,T)
      ENDIF
      XH=(H-H1)/(H2-H1)
      XS=(S-S1)/(S2-S1)
      NUHS9N=XH-XS
C
      END
C
C*******************************************************************
      SUBROUTINE WNHS3(XA,XB,F,P,T,EPS,X,IX)
C*******************************************************************
      IMPLICIT double precision(A-H,O-Z)
C
      EXTERNAL F
      X1=XA
      F1=F(X1,T,P)
      X3=XB
      F3=F(X3,T,P)
      IX=0
C++ SCHLEIFE ZUR ERMITTLUNG DER NULLSTELLE +++++++++++++++++++++++++++
      DO 100 I=1,200
      IF( F1 .NE. F3 ) THEN
       X=X1+(X3-X1)*F1/(F1-F3)
      ELSE
       IX=3
       GOTO 999
      END IF
      IF(X .LT. 0.D0) X=(X1+X3)/2.D0
      IF(DABS(X) .LT. 1.D-8) THEN
        IF(DABS(X-X1).LT.EPS) RETURN
      ELSE
        IF(DABS((X-X1)/X).LT.EPS) RETURN
      END IF
      F2=F(X,T,P)
      X2=X1-(X1-X3)/2.D0
      P1=F2*F1
      P3=F2*F3
      IF ((P1 .LT. 0.D0) .AND. (P3 .LT. 0.D0)) THEN
        P1=DABS(P1)
        P3=DABS(P3)
      ENDIF
      IF(P1 .LE. P3) THEN
        X3=X1
        F3=F1
      END IF
      X1=X
      F1=F2
      IF((X2-X3)*(X2-X1) .GE. 0.D0) GOTO 100
c      X=(X1+X3)/2.D0
      IF ((DABS(F1/F3)) .LE. (0.3D0)) THEN
         X=X1+DABS(F1/F3)*1.5D0*(X3-X1)
         IREM=1
      ELSE IF ((DABS(F1/F3)) .GE. (3.D0)) THEN
         IREM=2
      ELSE
         X=(X1+X3)/2.D0
         IREM=0
      ENDIF
      F2=F(X,T,P)
      IF ((((F2*F1) .GE. 0.D0) .AND. (IREM .EQ. 1)) .OR.
     1    (((F2*F3) .GE. 0.D0) .AND. (IREM .EQ. 2))) THEN
         X=(X1+X3)/2.D0
         F2=F(X,T,P)
      ENDIF
      IF(F2*F1 .LE. F2*F3) THEN
        X3=X1
        F3=F1
      END IF
      X1=X
      F1=F2
 100  CONTINUE
C++ ENDE DER SCHLEIFE ++++++++++++++++++++++++++++++++++++++++++++++++
      IX=1
 999  RETURN
      END
C***********************************************************************
      double precision FUNCTION PSATS2N(S)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      COMMON/CGR/PGR,TGR13
      COMMON/CCRIT/VC,SC,HC
      COMMON/CSUB2/TC,PC,DC
      double precision NUS2SATN
      EXTERNAL NUS2SATN
C
      EPS = 1.D-9
C
C     'ARBITRARY' BORDER FOR NEAR CRITICAL AREA
      IF (S .LE. 4.4123D0) THEN
          X = (S - SC) / (4.4123D0 - SC)
C         22.0639999990887 = SATURATION PRESSURE AT S'' = 4.4123
          PSATS2N = PC + X * (22.0639999990887D0 - PC)
          GOTO 999
      ENDIF
C
C     PSAT(350°C)=16.529 MPa
      PG = PSATTN(TGR13)
C     S(16.5;350°C)=5.21089kJ
      SG = SPT2N(PG,TGR13)                                
C     S''(21.5MPa)
      IF (S .LE. 4.71660853D+00) THEN                      
         P1 = 21.500D+00
         P2 = 22.064D+00
C     S''(21MPa)
      ELSEIF (S .LE. 4.80624D+00) THEN                    
         P1 = 21.0D+00
         P2 = 21.5D+00
      ELSEIF (S .LE. SG) THEN
         P1 = 18.1D+00                                                                
         P2 = 21.D0                                                              
C     S''(10MPa)
      ELSEIF (S .LE. 5.61589D+00) THEN                    
         P1 = 10.D+00
         P2 = 18.1D+00
C     S''(1MPa)
      ELSEIF (S .LE. 6.585D0) THEN                 
         P1 = 1.D0
         P2 = 10.D0
C     S''(0.1MPa)
      ELSEIF (S .LE. 7.359D0) THEN                      
         P1 = 0.1D0
         P2 = 1.D0
C     S''(0.01MPa)
      ELSEIF (S .LE. 8.149D0) THEN                       
         P1 = 0.01D0
         P2 = 0.1D0
      ELSE
         P1 = 0.0005D0
         P2 = 0.01D0
      ENDIF
      DUM = 1.D0
C
      CALL WNHSAT3(P1,P2,NUS2SATN,DUM,S,EPS,X,IX)
C
      PSATS2N = X
C
999   CONTINUE
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUS2SATN(P,S,DUM)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      Z=DUM
      IF (P .GT. 22.064D0) THEN
         PN=22.064D0
      ELSE   
         PN=P
      ENDIF
      T=TSATPN(PN)
      IF (T .LE. 623.15D0) THEN
         SN=SPT2N(PN,T)
      ELSE
         CALL FSATP(DV,DL,TOUT,T,PN)
         VN=1.D0/DV
         SN=SVT3N(VN,T)
      ENDIF
      NUS2SATN=SN-S
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION PSATH1N(H)
C***********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      double precision NUH1SATN
      EXTERNAL NUH1SATN
C
      EPS=1.D-9
C
      IF (H .LE. 417.5D0) THEN
         P1=0.001d0
         P2=0.1d0
      ELSEIF (H .LE. 1671.D0) THEN
         P1=0.1d0
         P2=10.d0
      ELSE
         P1=18.D0
         P2=21.D0
      ENDIF
      DUM=1.D0
C
      CALL WNHSAT3(P1,P2,NUH1SATN,DUM,H,EPS,X,IX)
C
      PSATH1N=X
C
      END      
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUH1SATN(P,H,DUM)
C***********************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      Z=DUM
      IF (P .GT. 22.064D0) THEN
         PN=22.064D0
      ELSE   
         PN=P
      ENDIF
      T=TSATPN(PN)
      IF (T .LE. 623.15D0) THEN
         HN=HPT1N(PN,T)
      ELSE
         CALL FSATP(DV,DL,TOUT,T,PN)
         VN=1.D0/DL
         HN=HVT3N(VN,T)
      ENDIF
      NUH1SATN=HN-H
      END
C
c******************************************************************
c******************************************************************
      subroutine fsatp(dvout,dlout,tout,tin,pin)
c******************************************************************
c
c     calculating tsat,saturated vapour and liquid volume for a given p
c
      implicit double precision(a-h,o-z)
c
      COMMON/CSUB2/TC,PC,DC
c
      p=pin
      DLOUT=0.d0
      DVOUT=0.d0
      TOUT=0.d0
c
      if( dabs(pc - pin) .lt. 1.d-05 ) then
          tout=tsatpn(p)
          DLOUT=dlest(tout)
          DVOUT=dvest(tout)
         return
      end if
c
      eps=1.d-09
c
      t=tin
      dle=dlest(t)
      dve=dvest(t)
c
         dlout=diter3(p,t,dle,eps)
       dvout=diter3(p,t,dve,eps)
        tout=t
c
      end
C
C***********************************************************************
      DOUBLEPRECISION FUNCTION FNR3N(V,T)
C***********************************************************************
C
C  First argument is specific volume in m^3/kg
C  Second argument is temperature in K
C  Returns the dimensionless reduced free energy
C
      IMPLICIT double precision (A-H,O-Z)
C
      PARAMETER (b01=  0.106 580 700 285 13 D01,
     1           b02= -0.157 328 452 902 39 D02,
     2           b03=  0.209 443 969 743 07 D02,
     3           b04= -0.768 677 078 787 16 D01,
     4           b05=  0.261 859 477 879 54 D01,
     5           b06= -0.280 807 811 486 20 D01,
     6           b07=  0.120 533 696 965 17 D01,
     7           b08= -0.845 668 128 125 02 D-02,
     8           b09= -0.126 543 154 777 14 D01,
     9           b10= -0.115 244 078 066 81 D01,
     1           b11=  0.885 210 439 843 18 D0,
     2           b12= -0.642 077 651 816 07 D0,
     3           b13=  0.384 934 601 866 71 D0,
     4           b14= -0.852 147 088 242 06 D0,
     5           b15=  0.489 722 815 418 77 D01,
     6           b16= -0.305 026 172 569 65 D01,
     7           b17=  0.394 205 368 791 54 D-01,
     8           b18=  0.125 584 084 243 08 D0,
     9           b19= -0.279 993 296 987 10 D0,
     1           b20=  0.138 997 995 694 60 D01,
     2           b21= -0.201 899 150 235 70 D01,
     3           b22= -0.821 476 371 739 63 D-02,
     4           b23= -0.475 960 357 349 23 D0,
     5           b24=  0.439 840 744 735 00 D-01,
     6           b25= -0.444 764 354 287 39 D0,
     7           b26=  0.905 720 707 197 33 D0,
     8           b27=  0.705 224 500 879 67 D0,
     9           b28=  0.107 705 126 263 32 D0,
     1           b29= -0.329 136 232 589 54 D0,
     2           b30= -0.508 710 620 411 58 D0,
     3           b31= -0.221 754 008 730 96 D-01,
     4           b32=  0.942 607 516 650 92 D-01,
     5           b33=  0.164 362 784 479 61 D0,
     6           b34= -0.135 033 722 413 48 D-01,
     7           b35= -0.148 343 453 524 72 D-01,
     8           b36=  0.579 229 536 280 84 D-03,
     9           b37=  0.323 089 047 037 11 D-02,
     1           b38=  0.809 648 029 962 15 D-04,
     2           b39= -0.165 576 797 950 37 D-03,
     3           b40= -0.449 238 990 618 15 D-04)
C
      PARAMETER ( tn= 647.096D0,
     2            rhonq= 1.D0/322.0D0)
C
      rho= 1.d0/v
      tau= tn/t 
      del= rho*rhonq
C
      tau2= tau*tau
      tau4= tau2*tau2
      tau12= tau4*tau4*tau4
C
      del2= del*del
C
      fnr3n= b01*dlog(del) + b02 + tau*(b03 + tau*(b04 +
     1     tau4*(tau*b05 + tau4*(b06 + tau2*b07 + tau12*(tau*b08 +
     2     del2*(b17 + tau4*(b18 + del*(b23 + del*(b27 +
     3     del*(b30 + del*(b33 + del2*(b35 + del*(b37 +
     4     del2*b40))) ))) ))) ))) ) +
     2     del* (tau2*(b09 + tau4*b10 + tau12*tau*(b11 + tau2*b12)) +
     4     del* (b13 + tau2*(b14 + tau4*(b15 + tau*b16)) +
     6     del* (b19 + tau2*(b20 + tau2*(b21 + tau12*b22)) +
     7     del* (b24 + tau2*(b25 + tau2*b26) +
     8     del* (tau*(b28 + tau2*b29) +
     1     del* (b31 + tau2*b32 +
     2     del* (tau2*b34 +
     3     del2*(tau2*b36 +
     4     del* (b38 + tau*b39))) ))) )))
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION DLest(T)
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      DATA TC/647.096d0/, DC/322.d0/,
     &     B1/1.99274064d0/, B2/1.09965342d0/,
     &     B3/-0.510839303d0/, B4/-1.75493479d0/,
     &     B5/-45.5170352d0/, B6/-6.74694450D+05/
C
      TAU = 1.D0 - T/TC
      TAUH = TAU**(1.D0/3.D0)
C
      DLest = DC* (1.D0 + B1*TAUH + B2*TAUH*TAUH + B3*TAUH**5.D0
     &        + B4*TAUH**16.D0 + B5*TAUH**43.D0 + B6*TAUH**110.D0)
C
      END
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION DVest(T)
C***********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      DATA TC/647.096d0/, DC/322.d0/,
     &     C1/-2.03150240d0/, C2/-2.68302940d0/,
     &     C3/-5.38626492d0/, C4/-17.2991605d0/,
     &     C5/-44.7586581d0/, C6/-63.9201063d0/
C
      TAU = 1.D0 - T/TC
      TAUH = TAU**(1.D0/6.D0)
C
      DH = C1*TAUH**2.D0 + C2*TAUH**4.D0 + C3*TAUH**8.D0
     &     + C4*TAUH**18.D0 + C5*TAUH**37.D0 + C6*TAUH**71.D0
      DVest = DC*DEXP(DH)
C
      END
C
C
      FUNCTION G0PT(P,T)
C
      IMPLICIT double precision (A-H,O-Z)
C
      PARAMETER (a01=-0.969 276 865 002 17 D01,
     1           a02= 0.100 866 559 680 18 D02,
     2           a03=-0.560 879 112 830 20 D-02,
     3           a04= 0.714 527 380 814 55 D-01,
     4           a05=-0.407 104 982 239 28 D0,
     5           a06= 0.142 408 191 714 44 D01,
     6           a07=-0.438 395 113 194 50 D01,
     7           a08=-0.284 086 324 607 72 D0,
     8           a09= 0.212 684 637 533 07 D-01)
C
      PARAMETER ( tn= 540.0D0,
     1            tnq= 1.D0/tn,
     2            r= 0.461526D0)
C
      tau= tn/t
      tauinv= t*tnq
C
      pi= p
C
      g0pt= (DLOG(pi) + a01 + tau*(a02 + tau*(a08 + tau*a09)) +
     1     tauinv*(a07 + tauinv*(a06 + tauinv*(a05 +
     2     tauinv*(a04 + tauinv*a03))))) * t*r
C
      END
C
C*************************************************************************************
      DOUBLE PRECISION FUNCTION dhdp1(P,T)
C*************************************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      parameter(a09=-0.16303507737913D+01,
     &          a10= 0.27182613947704D+01,
     2          a11= 0.12147472571259D+02,
     3          a13=-0.13971601220485D+02,
     4          a14=-0.10139813560979D+00,
     6          a15= 0.90547896843533D+00*2.D0,
     7          a17= 0.30487803863262D-01*2.D0,
     8          a18=-0.84709309503347D-02*2.D0,
     9          a19=-0.79051996435264D-11*2.D0,
     1          a20= 0.81058711826910D-01*3.D0,
     2          a22=-0.32702156038568D-05*3.D0,
     3          a23= 0.71724465059049D-02*4.D0,
     4          a24= 0.83376808703816D-03*4.D0,
     5          a25=-0.91740466143441D-09*4.D0,
     6          a26= 0.20734169140086D-02*5.D0,
     7          a27= 0.89603964175208D-05*8.D0,
     8          a28= 0.66877530790508D-06*8.D0,
     9          a29= 0.12755771455453D-13*21.D0,
     &          a30=-0.28710377452428D-15*23.D0,
     1          a31=-0.64016099916299D-18*29.D0,
     2          a32= 0.29806124175367D-18*30.D0,
     3          a33=-0.46640228230284D-19*31.D0,
     4          a34= 0.24531669269267D-20*32.D0)
c
        targ = 1386.d0/t-1.222d0
      targ2 = targ*targ
      targ4 = targ2*targ2
      targ8 = targ4*targ4
      targinv = 1.d0/targ
      targinv2 = targinv*targinv
      targinv4 = targinv2*targinv2
      targinv8 = targinv4*targinv4
      targinv9 = targinv8*targinv
C
      parg = 7.1d0-p*0.60496067755596D-01
      parg2 = parg*parg
      parg4 = parg2*parg2
      parg8 = parg4*parg4
C
      pt = parg*targinv
c
      hpt1z = (((((((((((pt*a34 + a33)*pt + a32)*pt + a31)*
     &   parg4*parg2*targinv9 + targinv2*a30)*parg2 +
     &   a29)*parg8*parg4*parg*targinv8*targinv8*targinv8 +
     &   targinv*a28)*targinv2 + targinv8*a27)*parg2*parg*
     &   targinv4 + targinv9*a26)*parg + targinv2*(targinv*
     &   a24 + targinv4*a23) + targ8*targ*a25)*parg + targ4*targ*
     &   a22 + targinv4*targinv*a20)*parg + targ8*targ8*a19 +
     &   targ2*a18 + targinv4*a15 + a17)*parg + targinv8*
     &   (targinv2*a09 + a10) + targ2*a14 + targinv2*a11 + a13
C
      dhdp1=hpt1z/(-16.53D0)
      END
C
C*************************************************************************************
      DOUBLE PRECISION FUNCTION dHdP2(P,T)
C*************************************************************************************
C
      IMPLICIT double precision (A-H,O-Z)
C
      parameter(a04=-0.44448764333452D+01,
     5          a05=-0.22926624714607D+02,
     6          a06=-0.43051902051180D+02,
     7          a07=-0.75253615672206D+02,
     8          a08=-0.82325284089205D-02*2.d0,
     9          a09=-0.94450864454513D-01*2.d0,
     &          a10=-0.39270508365637D+01*2.d0,
     1          a11=-0.76407372741773D+02*2.d0,
     2          a12=-0.23932578946761D+00*2.d0,
     4          a14= 0.10933624938123D-03*3.d0,
     5          a15=-0.24133119369638D-01*3.d0,
     6          a16=-0.22480892468696D+01*3.d0,
     7          a17=-0.35474272584198D+03*3.d0,
     8          a18=-0.19650645031516D-06*4.d0,
     9          a19= 0.63755087552931D-05*4.d0,
     &          a20= 0.36056766658235D-03*4.d0,
     1          a21= 0.39989127290421D-02*5.d0,
     2          a22=-0.12497164867770D-07*6.d0,
     3          a23=-0.84423037834823D+01*6.d0,
     4          a24=-0.20843876702652D+06*6.d0,
     6          a26=-0.34602240265361D-02*7.d0,
     7          a27=-0.24266223542696D+03*7.d0,
     8          a28= 0.22442547762780D-07*8.d0,
     9          a29=-0.73850273699100D+05*8.d0,
     &          a30= 0.64181736525088D-04*9.d0,
     1          a31= 0.10374663655276D-15*10.d0,
     2          a32=-0.25507450196258D-09*10.d0,
     3          a33=-0.34954795937691D-05*10.d0,
     4          a34=-0.58458099253864D-06*16.d0,
     5          a35= 0.13324803024175D+04*16.d0,
     6          a36=-0.47819819876448D+04*18.d0,
     7          a37= 0.44454513380586D-20*20.d0,
     8          a38= 0.26717467330171D-08*20.d0,
     9          a39=-0.50246518510642D-01*20.d0,
     &          a40=-0.30908182839692D-21*21.d0,
     1          a41= 0.49965138936998D-01*22.d0,
     2          a42=-0.12410752785137D-10*23.d0,
     3          a43= 0.47359492924764D-24*24.d0,
     4          a44= 0.55242716940682D-12*24.d0,
     5          a45=-0.13641135821518D-01*24.d0)
C
      tau= 540.d0/t
C      ti = t*0.18518518518518D-02
C
      targ = tau - 0.5d0
      targ2 = targ*targ
      targ3 = targ2*targ
      targ4 = targ3*targ
      targ6 = targ4*targ2
      targ7 = targ6*targ
      targ9 = targ7*targ2
      targ10 = targ9*targ
      targ13 = targ10*targ3
      targ14 = targ13*targ
      targ15 = targ14*targ
      targ19 = targ15*targ4
      targ20 = targ19*targ
      targ29 = targ20*targ9
C
      pi = p
      pi2 = pi*pi
      pi4 = pi2*pi2
C
      hpt2z=((((((((((((((((targ14*targ4*a45 + a44)*targ14 + a43)*
     &   pi2*targ6 + (pi*a42 + targ14*a41)*targ19)*pi +
     &   targ*a40)*pi + (targ13*a39 + a38)*targ15 + a37)*
     &   pi4 + (targ7*pi2*a36 + a35)*targ20*targ10 + targ9*a34)*
     &   pi4*pi2*targ19 + targ13*a33 + targ9*a32 +
     &   targ3*a31)*pi + targ10*targ2*a30)*pi + (targ14*targ14*a29 +
     &   a28)*targ7)*pi + (targ14*a27 + a26)*targ10)*pi + (targ19*
     &   a24 + a23)*targ15 + targ2*a22)*pi + targ6*
     &   a21)*pi + targ2*a20 + targ*a19 + a18)*pi + ((targ29*
     &   a17 + a16)*targ3 + a15)*targ2 + a14)*pi + (((targ29*
     &   a12 + a11)*targ3 + a10)*targ2 + a09)*targ + a08)*pi + (
     &   targ3*a07 + a06)*targ2 + targ*a05 + a04)
      dhdp2=hpt2z
c
      END
c
C*********************************************************************
      DOUBLE PRECISION FUNCTION dHdP5(P,T)
C*********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      PARAMETER(A7= -.63611419075956D+01,
     &          A8= .10049580368421D+01, 
     &          A10=-.16500637020131D-01*2.d0,
     &          A11=.17887679267013D-03*3.d0)
C
      TAU=T*1.D-3
      TAUM1=1.D0/TAU
      TAUM2=TAUM1*TAUM1
C
      dhdp5=((A10*TAUM2*TAUM2*TAUM2+A11*P)*P+A7)*TAUM2+A8
C
      END
c
c
C*********************************************************************
      DOUBLE PRECISION FUNCTION dHdV3(V,T)
C*********************************************************************
c
      IMPLICIT double precision (A-H,O-Z)
c
      parameter (f9 = -0.17520886815502D+01,
     &  f10 = -0.37231696861705D+01,
     &  f11 =  0.65367621353450D+01,
     &  f12 = -0.53340395459772D+01,
     &  f13 =  0.35531465412227D+00*2.d0,
     &  f14 = -0.15731521481920D+01*2.d0,
     &  f15 =  0.18081584968717D+02*2.d0,
     &  f16 = -0.12669975838924D+02*2.d0,
     &  f17 =  0.43664646488852D+00*2.d0,
     &  f18 =  0.16228889618024D+01*2.d0,
     &  f19 = -0.38767255915581D+00*3.d0,
     &  f20 =  0.32075594480472D+01*3.d0,
     &  f21 = -0.65227195048178D+01*3.d0,
     &  f22 = -0.72035213749267D-01*3.d0,
     &  f23 = -0.63703743166930D+01*3.d0,
     &  f24 =  0.81199175821825D-01*4.d0,
     &  f25 = -0.12316218802611D+01*4.d0,
     &  f26 =  0.33441092408796D+01*4.d0,
     &  f27 =  0.97643832897896D+01*4.d0,
     &  f28 =  0.29825229662283D+00*5.d0,
     &  f29 = -0.12152394310570D+01*5.d0,
     &  f30 = -0.72782785116784D+01*5.d0,
     &  f31 = -0.61407144380140D-01*6.d0,
     &  f32 =  0.34803030138386D+00*6.d0,
     &  f33 =  0.24274463510316D+01*6.d0,
     &  f34 = -0.56089416393544D-01*7.d0,
     &  f35 = -0.23277882648693D+00*8.d0,
     &  f36 =  0.29406244005770D-02*9.d0,
     &  f37 =  0.52189898432996D-01*9.d0,
     &  f38 =  0.37367361667631D-03*10.d0,
     &  f39 = -0.84059796975928D-03*10.d0,
     &  f40 = -0.76714125522093D-03*11.d0)
C
      dn=1.d0/(v*322.d0)
      dn2=dn*dn
c
      tn=647.096d0/t
      tn2=tn*tn
      tn4=tn2*tn2
      tn5=tn4*tn
      tn8=tn4*tn4
      tn12=tn8*tn4
C
      HVT3z =(((((((((((((dn2*f40 + f37)*dn + f35)*dn2 + f33)*dn +
     &   f30)*dn + f27)*dn + f23)*dn + f18)*tn4 + f17)*dn*tn5 +
     &   f12)*tn2 + f11)*tn*tn8 + f10)*tn4 + f9)*tn2 + ((((((((tn*f39
     &   + f38)*dn + tn2*f36)*dn2 + tn2*f34)*dn + tn2*f32 + f31)*
     &   dn + (tn2*f29 + f28)*tn)*dn + (tn2*f26 + f25)*
     &   tn2 + f24)*dn + ((tn12*f22 + f21)*tn2 + f20)*tn2 + f19)*
     &   dn + ((tn*f16 + f15)*tn4 + f14)*tn2 + f13)*dn)*t
         dhdv3=hvt3z/322.d0
c
      END
c
c
C*********************************************************************
      DOUBLE PRECISION FUNCTION dPdV3(V,T)
C*********************************************************************
c
      IMPLICIT double precision (A-H,O-Z)
c
      parameter (f1= 0.49189764279794D-03,
     &   f9 = -0.18137563991203D-05*2.d0,
     &  f10 = -0.16518055395610D-05*2.d0,
     &  f11 =  0.12687814703698D-05*2.d0,
     &  f12 = -0.92029667804989D-06*2.d0,
     &  f13 =  0.11034616587648D-05*3.d0,
     &  f14 = -0.24427828388075D-05*3.d0,
     &  f15 =  0.14038497646519D-04*3.d0,
     &  f16 = -0.87439446783466D-05*3.d0,
     &  f17 =  0.11300374350117D-06*3.d0,
     &  f18 =  0.36000198797746D-06*3.d0,
     &  f19 = -0.12039520470677D-05*4.d0,
     &  f20 =  0.59768188472930D-05*4.d0,
     &  f21 = -0.86815255166165D-05*4.d0,
     &  f22 = -0.35322922727656D-07*4.d0,
     &  f23 = -0.20465970175711D-05*4.d0,
     &  f24 =  0.25217135348393D-06*5.d0,
     &  f25 = -0.25499417810788D-05*5.d0,
     &  f26 =  0.51927162125460D-05*5.d0,
     &  f27 =  0.40432228943229D-05*5.d0,
     &  f28 =  0.77187447366156D-06*6.d0,
     &  f29 = -0.23587721876106D-05*6.d0,
     &  f30 = -0.36457015185726D-05*6.d0,
     &  f31 = -0.19070541732963D-06*7.d0,
     &  f32 =  0.81062958396863D-06*7.d0,
     &  f33 =  0.14134974870137D-05*7.d0,
     &  f34 = -0.13548168211001D-06*8.d0,
     &  f35 = -0.17009779063714D-06*9.d0,
     &  f36 =  0.74719422939563D-08*10.d0,
     &  f37 =  0.41677824835579D-07*10d0,
     &  f38 =  0.11604770704233D-08*11.d0,
     &  f39 = -0.23732297282870D-08*11.d0,
     &  f40 = -0.70828888764732D-09*12.d0)
C
      d=1.d0/v
      dn=d*0.31055900621118d-2
      dn2=dn*dn
c
      tn=647.096d0/t
      tn2=tn*tn
      tn4=tn2*tn2
      tn8=tn4*tn4
C
      dPdV3 = (((((((((((((dn2*f40 + f37)*dn + f35)*dn2 + f33)*dn +
     &   f30)*dn + f27)*dn + f23)*dn + f18)*tn4 + f17)*dn*tn*tn4 +
     &   f12)*tn2 + f11)*tn*tn8 + f10)*tn4 + f9)*tn2 + ((((((((tn*
     &   f39 + f38)*dn + tn2*f36)*dn2 + tn2*f34)*dn + tn2*f32 + f31)*
     &   dn + (tn2*f29 + f28)*tn)*dn + (tn2*f26 + f25)*
     &   tn2 + f24)*dn + ((tn4*tn8*f22 + f21)*tn2 + f20)*tn2 + f19)*
     &   dn + ((tn*f16 + f15)*tn4 + f14)*tn2 + f13)*dn + v*f1)*t*d
c
      END
c
C************************************************************************
      SUBROUTINE VISCVT(P,V,T,IREG,VISC)
C************************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      PARAMETER (H0 = 1.000000D0,
     1           H1 = 0.978197D0,
     2           H2 = 0.579829D0,
     3           H3 =-0.202354D0)
      PARAMETER (H00= 0.5132047D0,
     5           H10= 0.3205656D0,
     6           H40=-0.7782567D0,
     7           H50= 0.1885447D0,
     8           H01= 0.2151778D0,
     9           H11= 0.7317883D0,
     &           H21= 1.2410440D0,
     1           H31= 1.4767830D0,
     2           H02=-0.2818107D0,
     3           H12=-1.0707860D0,
     4           H22=-1.2631840D0,
     5           H03= 0.1778064D0,
     6           H13= 0.4605040D0,
     7           H23= 0.2340379D0,
     8           H33=-0.4924179D0,
     9           H04=-0.0417661D0,
     &           H34= 0.1600435D0,
     1           H15=-0.01578386D0,
     2           H36=-0.003629481D0)
C
      IF (IREG .EQ. 0) THEN
         CALL REGSOPT(P,T,IREG)
      ENDIF
C
      IF (((IREG .NE. 0) .AND. (IREG .NE. 5) .AND. (IREG .NE. 9)) .OR.
     1((IREG .EQ. 5) .AND. (T .LE. 1173.15D0)))THEN
         D=1.D0/V
      ELSE 
c         WRITE(IANORM,*)'DATAPOINT OUTSIDE THE TEMPERATURE AND '
c         WRITE(IANORM,*)'PRESSURE RANGE OF THE VISCOSITY FUNCTION'
         VISC=0.D0
         GOTO 100
      ENDIF
C
      TR=T/6.47226D2
      DR=D/3.17763D2
      TRINV=1.D0/TR
      TRINV2=TRINV*TRINV
      
      TH=TRINV-1.D0
      TH2=TH*TH
      TH3=TH2*TH
C      
      DH=DR-1.D0
      DH2=DH*DH
      DH3=DH2*DH
C
      VISC0=DSQRT(TR)/(H0+(H1*TRINV)+(H2*TRINV2)+(H3*TRINV*TRINV2))
      VISC1=DEXP(DR*(H00+DH*(H01+TH*(H11+TH*(H21+TH*H31)))
     #           +TH*(H10+TH3*(H40+TH*H50))+DH2*(H02+TH*(H12+TH*H22))
     #           +DH3*(H03+TH*(H13+TH*(H23+TH*H33))
     #           +DH*(H04+H34*TH3+DH*(H15*TH+DH*H36*TH3)))))  
C
      VISC=VISC0*VISC1*55.071D-6
C
  100 CONTINUE
C      
      END
C
C************************************************************************
      SUBROUTINE THCONVT(P,V,T,IREG,THCON)
C************************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      PARAMETER (a0 = 1.000000D0,
     1           a1 = 6.978267D0,
     2           a2 = 2.599096D0,
     3           a3 =-0.998254D0)
      PARAMETER(B00 = 1.3293046D0,
     1          B01 =-0.40452437D0,
     2          B02 = 0.24409490D0,
     3          B03 = 0.018660751D0,
     4          B04 =-0.12961068D0,
     5          B05 = 0.044809953D0,
     6          B10 = 1.7018363D0,
     7          B11 =-2.2156845D0,
     8          B12 = 1.6511057D0,
     9          B13 =-0.76736002D0,
     &          B14 = 0.37283344D0,
     1          B15 =-0.11203160D0,
     2          B20 = 5.2246158D0,
     3          B21 =-10.124111D0,
     4          B22 = 4.9874687D0,
     5          B23 =-0.27297694D0,
     6          B24 =-0.43083393D0,
     7          B25 = 0.13333849D0,
     8          B30 = 8.7127675D0,
     9          B31 =-9.5000611D0,
     &          B32 = 4.3786606D0,
     1          B33 =-0.91783782D0,
     2          B40 =-1.8525999D0,
     3          B41 = 0.93404690D0)
C
      IF (IREG .EQ. 0) THEN
         CALL REGSOPT(P,T,IREG)
      ENDIF
C
      IF ((IREG .NE. 0) .AND. (IREG .NE. 5) .AND. (IREG .NE. 9)) THEN
         D=1.D0/V
      ELSE
         THCON=0.D0
         GOTO 100
      ENDIF
C
      TR=T/6.47226D2
      TR2=TR*TR
      TH=(1.D0/TR-1.D0)
      TH2=TH*TH
      TH3=TH2*TH
      TZ=TR-1.D0
      DR=D/3.17763D2
      DH=DR-1.D0
      DH2=DH*DH
      DH3=DH2*DH
C
      CALL CONVAL(T,D,P,IREG,VISC0,VISC1,DPDT,CHI)
C
      THCON0=DSQRT(TR)/(a0+a1/TR+a2/TR2+a3/TR/TR2)
      THCON1=DEXP(DR*(B00+DH*B01+DH2*B02+DH3*B03
     &            +DH2*DH2*B04
     &            +DH2*DH3*B05
     1        +TH*B10+TH*DH*B11+TH*DH2*B12+TH*DH3*B13
     &            +TH*DH2*DH2*B14
     &            +TH*DH2*DH3*B15
     1        +TH2*B20+TH2*DH*B21+TH2*DH2*B22+TH2*DH3*B23
     &            +TH2*DH2*DH2*B24
     &            +TH2*DH2*DH3*B25
     1        +TH3*B30+TH3*DH*B31+TH3*DH2*B32
     &            +TH3*DH3*B33
     1        +TH2*TH2*B40+TH2*TH2*DH*B41))
      THCON2=0.0013848D0/VISC0/VISC1*TR*TR/DR/DR*DPDT*DPDT
     1       *CHI**0.4678D0*DSQRT(DR)*DEXP((-18.66D0)*TZ*TZ-DH2*DH2)
C
      THCON=(THCON0*THCON1+THCON2)*0.4945D0
C
  100 CONTINUE
C      
      END
C
C************************************************************************
      SUBROUTINE CONVAL(T,D,P,IREG,VISC0,VISC1,DPDT,CHI)
C************************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      PARAMETER (H0 = 1.000000D0,
     1           H1 = 0.978197D0,
     2           H2 = 0.579829D0,
     3           H3 =-0.202354D0)
      PARAMETER (H00= 0.5132047D0,
     5           H10= 0.3205656D0,
     6           H40=-0.7782567D0,
     7           H50= 0.1885447D0,
     8           H01= 0.2151778D0,
     9           H11= 0.7317883D0,
     &           H21= 1.2410440D0,
     1           H31= 1.4767830D0,
     2           H02=-0.2818107D0,
     3           H12=-1.0707860D0,
     4           H22=-1.2631840D0,
     5           H03= 0.1778064D0,
     6           H13= 0.4605040D0,
     7           H23= 0.2340379D0,
     8           H33=-0.4924179D0,
     9           H04=-0.0417661D0,
     &           H34= 0.1600435D0,
     1           H15=-0.01578386D0,
     2           H36=-0.003629481D0)
C
      PARAMETER (c15= -0.471 843 210 732 67 D-03*2.D0,
     6           c16= -0.300 017 807 930 26 D-03*2.D0,
     7           c17=  0.476 613 939 069 87 D-04*2.D0,
     8           c18= -0.441 418 453 308 46 D-05*2.D0,
     9           c19= -0.726 949 962 975 94 D-15*2.D0,
     1           c20= -0.316 796 448 450 54 D-04*6.D0,
     2           c21= -0.282 707 979 853 12 D-05*6.D0,
     3           c22= -0.852 051 281 201 03 D-09*6.D0,
     4           c23= -0.224 252 819 080 00 D-05*12.D0,
     5           c24= -0.651 712 228 956 01 D-06*12.D0,
     6           c25= -0.143 417 299 379 24 D-12*12.D0,
     7           c26= -0.405 169 968 601 17 D-06*20.D0,
     8           c27= -0.127 343 017 416 41 D-08*56.D0,
     9           c28= -0.174 248 712 306 34 D-09*56.D0,
     1           c29= -0.687 621 312 955 31 D-18*420.D0,
     2           c30=  0.144 783 078 285 21 D-19*506.D0,
     3           c31=  0.263 357 816 627 95 D-22*812.D0,
     4           c32= -0.119 476 226 400 71 D-22*870.D0,
     5           c33=  0.182 280 945 814 04 D-23*930.D0,
     6           c34= -0.935 370 872 924 58 D-25*992.D0)
C
      PARAMETER (b09=  0.283 190 801 238 04 D-03,
     9           b10= -0.607 063 015 658 74 D-03,
     1           b11= -0.189 900 682 184 19 D-01,
     2           b12= -0.325 297 487 705 05 D-01,
     3           b13= -0.218 417 171 754 14 D-01,
     4           b14= -0.528 383 579 699 30 D-04,
     5           b15= -0.471 843 210 732 67 D-03*2.D0,
     6           b16= -0.300 017 807 930 26 D-03*2.D0,
     7           b17=  0.476 613 939 069 87 D-04*2.D0,
     8           b18= -0.441 418 453 308 46 D-05*2.D0,
     9           b19= -0.726 949 962 975 94 D-15*2.D0,
     1           b20= -0.316 796 448 450 54 D-04*3.D0,
     2           b21= -0.282 707 979 853 12 D-05*3.D0,
     3           b22= -0.852 051 281 201 03 D-09*3.D0,
     4           b23= -0.224 252 819 080 00 D-05*4.D0,
     5           b24= -0.651 712 228 956 01 D-06*4.D0,
     6           b25= -0.143 417 299 379 24 D-12*4.D0,
     7           b26= -0.405 169 968 601 17 D-06*5.D0,
     8           b27= -0.127 343 017 416 41 D-08*8.D0,
     9           b28= -0.174 248 712 306 34 D-09*8.D0,
     1           b29= -0.687 621 312 955 31 D-18*21.D0,
     2           b30=  0.144 783 078 285 21 D-19*23.D0,
     3           b31=  0.263 357 816 627 95 D-22*29.D0,
     4           b32= -0.119 476 226 400 71 D-22*30.D0,
     5           b33=  0.182 280 945 814 04 D-23*31.D0,
     6           b34= -0.935 370 872 924 58 D-25*32.D0)
C
      PARAMETER (f09= -0.283 190 801 238 04 D-03*9.D0,
     9           f10=  0.607 063 015 658 74 D-03*7.D0,
     1           f11=  0.189 900 682 184 19 D-01*1.D0,
     3           f13= -0.218 417 171 754 14 D-01,
     4           f14= -0.528 383 579 699 30 D-04*3.D0,
     5           f15=  0.471 843 210 732 67 D-03*6.D0,
     7           f17=  0.476 613 939 069 87 D-04*2.D0,
     8           f18= -0.441 418 453 308 46 D-05*6.D0,
     9           f19= -0.726 949 962 975 94 D-15*34.D0,
     1           f20=  0.316 796 448 450 54 D-04*12.D0,
     3           f22= -0.852 051 281 201 03 D-09*18.D0,
     4           f23=  0.224 252 819 080 00 D-05*20.D0,
     5           f24=  0.651 712 228 956 01 D-06*8.D0,
     6           f25= -0.143 417 299 379 24 D-12*40.D0,
     7           f26=  0.405 169 968 601 17 D-06*40.D0,
     8           f27=  0.127 343 017 416 41 D-08*88.D0,
     9           f28=  0.174 248 712 306 34 D-09*48.D0,
     1           f29=  0.687 621 312 955 31 D-18*609.D0,
     2           f30= -0.144 783 078 285 21 D-19*713.D0,
     3           f31= -0.263 357 816 627 95 D-22*1102.D0,
     4           f32=  0.119 476 226 400 71 D-22*1170.D0,
     5           f33= -0.182 280 945 814 04 D-23*1240.D0,
     6           f34=  0.935 370 872 924 58 D-25*1312.D0)
C
      PARAMETER (a06=-0.330 326 416 702 03 D-04*2.D0,
     6           a07=-0.189 489 875 163 15 D-03*2.D0,
     7           a08=-0.393 927 772 433 55 D-02*2.D0,
     8           a09=-0.437 972 956 505 73 D-01*2.D0,
     9           a10=-0.266 745 479 140 87 D-04*2.D0,
     1           a11= 0.204 817 376 923 09 D-07*6.D0,
     2           a12= 0.438 706 672 844 35 D-06*6.D0,
     3           a13=-0.322 776 772 385 70 D-04*6.D0,
     4           a14=-0.150 339 245 421 48 D-02*6.D0,
     5           a15=-0.406 682 535 626 49 D-01*6.D0,
     6           a16=-0.788 473 095 593 67 D-09*12.D0,
     7           a17= 0.127 907 178 522 85 D-07*12.D0,
     8           a18= 0.482 253 727 185 07 D-06*12.D0,
     9           a19= 0.229 220 763 376 61 D-05*20.D0,
     1           a20=-0.167 147 664 510 61 D-10*30.D0,
     2           a21=-0.211 714 723 213 55 D-02*30.D0,
     3           a22=-0.238 957 419 341 04 D02*30.D0,
     4           a23=-0.590 595 643 242 70 D-17*42.D0,
     5           a24=-0.126 218 088 991 01 D-05*42.D0,
     6           a25=-0.389 468 424 357 39 D-01*42.D0,
     7           a26= 0.112 562 113 604 59 D-10*56.D0,
     8           a27=-0.823 113 408 979 98 D01*56.D0,
     9           a28= 0.198 097 128 020 88 D-07*72.D0,
     1           a29= 0.104 069 652 101 74 D-18*90.D0,
     2           a30=-0.102 347 470 959 29 D-12*90.D0,
     3           a31=-0.100 181 793 795 11 D-08*90.D0,
     4           a32=-0.808 829 086 469 85 D-10*240.D0,
     5           a33= 0.106 930 318 794 09 D0*240.D0,
     6           a34=-0.336 622 505 741 71 D0*306.D0,
     7           a35= 0.891 858 453 554 21 D-24*380.D0,
     8           a36= 0.306 293 168 762 32 D-12*380.D0,
     9           a37=-0.420 024 676 982 08 D-05*380.D0,
     1           a38=-0.590 560 296 856 39 D-25*420.D0,
     2           a39= 0.378 269 476 134 57 D-05*462.D0,
     3           a40=-0.127 686 089 346 81 D-14*506.D0,
     4           a41= 0.730 876 105 950 61 D-28*552.D0,
     5           a42= 0.554 147 153 507 78 D-16*552.D0,
     6           a43=-0.943 697 072 412 10 D-06*552.D0)
C
      PARAMETER (e02=-0.178 348 622 923 58 D-01,
     2           e03=-0.459 960 136 963 65 D-01*2.D0,
     3           e04=-0.575 812 590 834 32 D-01*3.D0,
     4           e05=-0.503 252 787 279 30 D-01*6.D0,
     5           e06=-0.330 326 416 702 03 D-04*2.D0,
     6           e07=-0.189 489 875 163 15 D-03*4.D0,
     7           e08=-0.393 927 772 433 55 D-02*8.D0,
     8           e09=-0.437 972 956 505 73 D-01*14.D0,
     9           e10=-0.266 745 479 140 87 D-04*72.D0,
     2           e12= 0.438 706 672 844 35 D-06*3.D0,
     3           e13=-0.322 776 772 385 70 D-04*9.D0,
     4           e14=-0.150 339 245 421 48 D-02*18.D0,
     5           e15=-0.406 682 535 626 49 D-01*105.D0,
     6           e16=-0.788 473 095 593 67 D-09*4.D0,
     7           e17= 0.127 907 178 522 85 D-07*8.D0,
     8           e18= 0.482 253 727 185 07 D-06*12.D0,
     9           e19= 0.229 220 763 376 61 D-05*35.D0,
     1           e20=-0.167 147 664 510 61 D-10*18.D0,
     2           e21=-0.211 714 723 213 55 D-02*96.D0,
     3           e22=-0.238 957 419 341 04 D02*210.D0,
     5           e24=-0.126 218 088 991 01 D-05*77.D0,
     6           e25=-0.389 468 424 357 39 D-01*175.D0,
     7           e26= 0.112 562 113 604 59 D-10*64.D0,
     8           e27=-0.823 113 408 979 98 D01*288.D0,
     9           e28= 0.198 097 128 020 88 D-07*117.D0,
     1           e29= 0.104 069 652 101 74 D-18*40.D0,
     2           e30=-0.102 347 470 959 29 D-12*100.D0,
     3           e31=-0.100 181 793 795 11 D-08*140.D0,
     4           e32=-0.808 829 086 469 85 D-10*464.D0,
     5           e33= 0.106 930 318 794 09 D0*800.D0,
     6           e34=-0.336 622 505 741 71 D0*1026.D0,
     7           e35= 0.891 858 453 554 21 D-24*400.D0,
     8           e36= 0.306 293 168 762 32 D-12*700.D0,
     9           e37=-0.420 024 676 982 08 D-05*960.D0,
     1           e38=-0.590 560 296 856 39 D-25*441.D0,
     2           e39= 0.378 269 476 134 57 D-05*1166.D0,
     3           e40=-0.127 686 089 346 81 D-14*897.D0,
     4           e41= 0.730 876 105 950 61 D-28*624.D0,
     5           e42= 0.554 147 153 507 78 D-16*960.D0,
     6           e43=-0.943 697 072 412 10 D-06*1392.D0)
C
      PARAMETER (d01=-0.177 317 424 732 13 D-02,
     1           d02=-0.178 348 622 923 58 D-01,
     2           d03=-0.459 960 136 963 65 D-01,
     3           d04=-0.575 812 590 834 32 D-01,
     4           d05=-0.503 252 787 279 30 D-01,
     5           d06=-0.330 326 416 702 03 D-04*2.D0,
     6           d07=-0.189 489 875 163 15 D-03*2.D0,
     7           d08=-0.393 927 772 433 55 D-02*2.D0,
     8           d09=-0.437 972 956 505 73 D-01*2.D0,
     9           d10=-0.266 745 479 140 87 D-04*2.D0,
     1           d11= 0.204 817 376 923 09 D-07*3.D0,
     2           d12= 0.438 706 672 844 35 D-06*3.D0,
     3           d13=-0.322 776 772 385 70 D-04*3.D0,
     4           d14=-0.150 339 245 421 48 D-02*3.D0,
     5           d15=-0.406 682 535 626 49 D-01*3.D0,
     6           d16=-0.788 473 095 593 67 D-09*4.D0,
     7           d17= 0.127 907 178 522 85 D-07*4.D0,
     8           d18= 0.482 253 727 185 07 D-06*4.D0,
     9           d19= 0.229 220 763 376 61 D-05*5.D0,
     1           d20=-0.167 147 664 510 61 D-10*6.D0,
     2           d21=-0.211 714 723 213 55 D-02*6.D0,
     3           d22=-0.238 957 419 341 04 D02*6.D0,
     4           d23=-0.590 595 643 242 70 D-17*7.D0,
     5           d24=-0.126 218 088 991 01 D-05*7.D0,
     6           d25=-0.389 468 424 357 39 D-01*7.D0,
     7           d26= 0.112 562 113 604 59 D-10*8.D0,
     8           d27=-0.823 113 408 979 98 D01*8.D0,
     9           d28= 0.198 097 128 020 88 D-07*9.D0,
     1           d29= 0.104 069 652 101 74 D-18*10.D0,
     2           d30=-0.102 347 470 959 29 D-12*10.D0,
     3           d31=-0.100 181 793 795 11 D-08*10.D0,
     4           d32=-0.808 829 086 469 85 D-10*16.D0,
     5           d33= 0.106 930 318 794 09 D0*16.D0,
     6           d34=-0.336 622 505 741 71 D0*18.D0,
     7           d35= 0.891 858 453 554 21 D-24*20.D0,
     8           d36= 0.306 293 168 762 32 D-12*20.D0,
     9           d37=-0.420 024 676 982 08 D-05*20.D0,
     1           d38=-0.590 560 296 856 39 D-25*21.D0,
     2           d39= 0.378 269 476 134 57 D-05*22.D0,
     3           d40=-0.127 686 089 346 81 D-14*23.D0,
     4           d41= 0.730 876 105 950 61 D-28*24.D0,
     5           d42= 0.554 147 153 507 78 D-16*24.D0,
     6           d43=-0.943 697 072 412 10 D-06*24.D0)
C
      PARAMETER (r01=  0.106 580 700 285 13 D01,
     1           r09= -0.126 543 154 777 14 D01,
     9           r10= -0.115 244 078 066 81 D01,
     1           r11=  0.885 210 439 843 18 D0,
     2           r12= -0.642 077 651 816 07 D0,
     3           r13=  0.384 934 601 866 71 D0*2.D0,
     4           r14= -0.852 147 088 242 06 D0*2.D0,
     5           r15=  0.489 722 815 418 77 D01*2.D0,
     6           r16= -0.305 026 172 569 65 D01*2.D0,
     7           r17=  0.394 205 368 791 54 D-01*2.D0,
     8           r18=  0.125 584 084 243 08 D0*2.D0,
     9           r19= -0.279 993 296 987 10 D0*3.D0,
     1           r20=  0.138 997 995 694 60 D01*3.D0,
     2           r21= -0.201 899 150 235 70 D01*3.D0,
     3           r22= -0.821 476 371 739 63 D-02*3.D0,
     4           r23= -0.475 960 357 349 23 D0*3.D0,
     5           r24=  0.439 840 744 735 00 D-01*4.D0,
     6           r25= -0.444 764 354 287 39 D0*4.D0,
     7           r26=  0.905 720 707 197 33 D0*4.D0,
     8           r27=  0.705 224 500 879 67 D0*4.D0,
     9           r28=  0.107 705 126 263 32 D0*5.D0,
     1           r29= -0.329 136 232 589 54 D0*5.D0,
     2           r30= -0.508 710 620 411 58 D0*5.D0,
     3           r31= -0.221 754 008 730 96 D-01*6.D0,
     4           r32=  0.942 607 516 650 92 D-01*6.D0,
     5           r33=  0.164 362 784 479 61 D0*6.D0,
     6           r34= -0.135 033 722 413 48 D-01*7.D0,
     7           r35= -0.148 343 453 524 72 D-01*8.D0,
     8           r36=  0.579 229 536 280 84 D-03*9.D0,
     9           r37=  0.323 089 047 037 11 D-02*9.D0,
     1           r38=  0.809 648 029 962 15 D-04*10.D0,
     2           r39= -0.165 576 797 950 37 D-03*10.D0,
     3           r40= -0.449 238 990 618 15 D-04*11.D0)
      PARAMETER (s01=  0.106 580 700 285 13 D01,
     3           s13=  0.384 934 601 866 71 D0*2.D0,
     4           s14= -0.852 147 088 242 06 D0*2.D0,
     5           s15=  0.489 722 815 418 77 D01*2.D0,
     6           s16= -0.305 026 172 569 65 D01*2.D0,
     7           s17=  0.394 205 368 791 54 D-01*2.D0,
     8           s18=  0.125 584 084 243 08 D0*2.D0,
     9           s19= -0.279 993 296 987 10 D0*6.D0,
     1           s20=  0.138 997 995 694 60 D01*6.D0,
     2           s21= -0.201 899 150 235 70 D01*6.D0,
     3           s22= -0.821 476 371 739 63 D-02*6.D0,
     4           s23= -0.475 960 357 349 23 D0*6.D0,
     5           s24=  0.439 840 744 735 00 D-01*12.D0,
     6           s25= -0.444 764 354 287 39 D0*12.D0,
     7           s26=  0.905 720 707 197 33 D0*12.D0,
     8           s27=  0.705 224 500 879 67 D0*12.D0,
     9           s28=  0.107 705 126 263 32 D0*20.D0,
     1           s29= -0.329 136 232 589 54 D0*20.D0,
     2           s30= -0.508 710 620 411 58 D0*20.D0,
     3           s31= -0.221 754 008 730 96 D-01*30.D0,
     4           s32=  0.942 607 516 650 92 D-01*30.D0,
     5           s33=  0.164 362 784 479 61 D0*30.D0,
     6           s34= -0.135 033 722 413 48 D-01*42.D0,
     7           s35= -0.148 343 453 524 72 D-01*56.D0,
     8           s36=  0.579 229 536 280 84 D-03*72.D0,
     9           s37=  0.323 089 047 037 11 D-02*72.D0,
     1           s38=  0.809 648 029 962 15 D-04*90.D0,
     2           s39= -0.165 576 797 950 37 D-03*90.D0,
     3           s40= -0.449 238 990 618 15 D-04*110.D0)
      PARAMETER (t09= -0.126 543 154 777 14 D01*2.D0,
     9           t10= -0.115 244 078 066 81 D01*6.D0,
     1           t11=  0.885 210 439 843 18 D0*15.D0,
     2           t12= -0.642 077 651 816 07 D0*17.D0,
     4           t14= -0.852 147 088 242 06 D0*4.D0,
     5           t15=  0.489 722 815 418 77 D01*12.D0,
     6           t16= -0.305 026 172 569 65 D01*14.D0,
     7           t17=  0.394 205 368 791 54 D-01*44.D0,
     8           t18=  0.125 584 084 243 08 D0*52.D0,
     1           t20=  0.138 997 995 694 60 D01*6.D0,
     2           t21= -0.201 899 150 235 70 D01*12.D0,
     3           t22= -0.821 476 371 739 63 D-02*48.D0,
     4           t23= -0.475 960 357 349 23 D0*78.D0,
     6           t25= -0.444 764 354 287 39 D0*8.D0,
     7           t26=  0.905 720 707 197 33 D0*16.D0,
     8           t27=  0.705 224 500 879 67 D0*104.D0,
     9           t28=  0.107 705 126 263 32 D0*5.D0,
     1           t29= -0.329 136 232 589 54 D0*15.D0,
     2           t30= -0.508 710 620 411 58 D0*130.D0,
     4           t32=  0.942 607 516 650 92 D-01*12.D0,
     5           t33=  0.164 362 784 479 61 D0*156.D0,
     6           t34= -0.135 033 722 413 48 D-01*14.D0,
     7           t35= -0.148 343 453 524 72 D-01*208.D0,
     8           t36=  0.579 229 536 280 84 D-03*18.D0,
     9           t37=  0.323 089 047 037 11 D-02*234.D0,
     2           t39= -0.165 576 797 950 37 D-03*10.D0,
     3           t40= -0.449 238 990 618 15 D-04*286.D0)
C
      PARAMETER (R=0.461526D3)
C
      PARAMETER (TNR=647.226D0,
     1           DNR=317.763D0,
     2           PNR=22.115D0)
C
      TR=T/TNR
      DR=D/DNR
C
      TH=1.D0/TR-1.D0
      TH2=TH*TH
      TH3=TH2*TH
C      
      DH=DR-1.D0
      DH2=DH*DH
      DH3=DH2*DH
C
      VISC0=DSQRT(TR)/(H0+(H1/TR)+(H2/TR/TR)+(H3/TR/TR/TR))
      VISC1=DEXP(DR*(H00+H10*TH+H40*TH2*TH2+H50*TH2*TH3+
     1               H01*DH+H11*TH*DH+H21*TH2*DH+H31*TH3*DH+
     2               H02*DH2+H12*TH*DH2+H22*TH2*DH2+
     3               H03*DH3+H13*TH*DH3+H23*TH2*DH3+
     4               H33*TH3*DH3+H04*DH2*DH2+
     5               H34*TH3*DH2*DH2+H15*TH*DH2*DH3+
     6               H36*TH3*DH3*DH3))
C
      IF (IREG .EQ. 1) THEN
C
         TN= 1386.D0
         PN= 16.53D0
C
         tau= tn/t
         TAU2=TAU*TAU
         taun= tau - 1.222d0
         taun2= taun*taun
         taun3= taun2*taun
         taun6= taun3*taun3
         taun12= taun6*taun6
         tauninv= 1.d0/taun
         tauninv2= tauninv*tauninv
         tauninv3= tauninv2*tauninv
         tauninv7= tauninv3*tauninv3*tauninv
C
         PI=P/PN
         pin= 7.1d0 - p/pn
         pin2= pin*pin
         pin3= pin2*pin
         pin6= pin3*pin3
C         
         gapp= tauninv3*(c15 + taun3*(c16 + taun*(c17 +
     5       taun2*(c18 + taun12*taun2*c19))) +
     6       pin*  tauninv*(c20 + taun3*taun*(c21 + taun6*c22) +
     7       pin*  tauninv*(c23 + taun3*(c24 + taun12*c25) +
     8       pin*  tauninv3*(c26 +
     9       pin3* (tauninv3*(c27 + taun3*taun2*c28) +
     1       pin6*pin6*pin*tauninv7*tauninv7*tauninv7*(c29 +
     2       pin2* tauninv2*(c30 +
     3       pin6* tauninv7*(c31 +
     4       pin*  tauninv*(c32 +
     5       pin*  tauninv*(c33 +
     6       pin*  tauninv*c34))) ))) ))) )
C
         gap= -(tauninv7*(tauninv2*b09 + b10 + taun6*(b11 +
     1      taun*(b12 + taun*(b13 + taun2*b14)))) +
     2      pin*  tauninv3*(b15 + taun3*(b16 + taun*(b17 +
     3            taun2*(b18 + taun12*taun2*b19))) +
     4      pin*  tauninv*(b20 + taun3*taun*(b21 + taun6*b22) +
     5      pin*  tauninv*(b23 + taun3*(b24 + taun12*b25) +
     6      pin*  tauninv3*(b26 +
     7      pin3* (tauninv3*(b27 + taun3*taun2*b28) +
     8      pin6*pin6*pin*tauninv7*tauninv7*tauninv7*(b29 +
     9      pin2* tauninv2*(b30 +
     1      pin6* tauninv7*(b31 +
     2      pin*  tauninv*(b32 +
     3      pin*  tauninv*(b33 +
     4      pin*  tauninv*b34))) ))) ))) ))
C
         gapt= -(tauninv*(tauninv7*(tauninv2*f09 + f10 + taun6*(f11 +
     3       taun2*(f13 + taun2*f14))) +
     4       pin*  tauninv3*(f15 + taun3*taun*(f17 +
     5             taun2*(f18 + taun12*taun2*f19)) +
     6       pin*  tauninv*(f20 + taun6*taun3*taun*f22 +
     7       pin*  tauninv*(f23 + taun3*(f24 + taun12*f25) +
     8       pin*  tauninv3*(f26 +
     9       pin3*(tauninv3*(f27 + taun3*taun2*f28) +
     1       pin6*pin6*pin*tauninv7*tauninv7*tauninv7*(f29 +
     2       pin2* tauninv2*(f30 +
     3       pin6* tauninv7*(f31 +
     4       pin*  tauninv*(f32 +
     5       pin*  tauninv*(f33 +
     6       pin*  tauninv*f34))) ))) ))) )))
C     
         GA0=0.D0
C
      ELSE IF (IREG .EQ. 2) THEN
C
         TN= 540.D0
         PN= 1.D0
C
         pi= p/pn
         tau=tn/t
         TAU2=TAU*TAU
C
         taun= tau - 0.5d0
         taun2= taun*taun
         taun3= taun2*taun
         taun4= taun3*taun
         taun6= taun4*taun2
         taun7= taun6*taun
         taun9= taun7*taun2
         taun10= taun9*taun
         taun13= taun10*taun3
         taun14= taun13*taun
         taun15= taun14*taun
         taun19= taun15*taun4
         taun20= taun19*taun
         taun29= taun20*taun9
C
         pi2= pi*pi
         pi4= pi2*pi2
C
         gap= d01 + taun*(d02 + taun*(d03 + taun*(d04 + taun3*d05))) +
     3      pi* (taun*(d06 + taun*(d07 + taun2*(d08 +
     4         taun3*(d09 + taun29*d10)))) +
     5      pi* (d11 + taun*(d12 + taun2*(d13 + taun3*(d14 +
     6         taun29*d15))) +
     7      pi* (taun*(d16 + taun*(d17 + taun*d18)) +
     8      pi* (taun7*d19 +
     9      pi* (taun3*(d20 + taun13*(d21 + taun19*d22)) +
     1      pi* (d23 + taun10*taun*(d24 + taun14*d25) +
     2      pi* (taun7*(taun*d26 + taun29*d27) +
     3      pi* (taun13*d28 +
     4      pi* (taun4*(d29 + taun6*(d30 + taun4*d31)) +
     5      pi4*pi2*taun20*(taun9*d32 + taun20*taun10*(d33 +
     6         pi2*taun7*d34) +
     7      pi4*(d35 + taun15*(d36 + taun13*d37) +
     8      pi* (taun*d38 +
     9      pi* (taun19*(taun14*d39 +
     1          pi* d40) +
     2      pi2*(taun6*(d41 + taun14*(d42 +
     3      taun14*taun4*d43))) ))) ))) ))) ))) )
C
         gapp= (taun*(a06 + taun*(a07 + taun2*(a08 +
     4      taun3*(a09 + taun29*a10)))) +
     5      pi* (a11 + taun*(a12 + taun2*(a13 + taun3*(a14 +
     6         taun29*a15))) +
     7      pi* (taun*(a16 + taun*(a17 + taun*a18)) +
     8      pi* (taun7*a19 +
     9      pi* (taun3*(a20 + taun13*(a21 + taun19*a22)) +
     1      pi* (a23 + taun10*taun*(a24 + taun14*a25) +
     2      pi* (taun7*(taun*a26 + taun29*a27) +
     3      pi* (taun13*a28 +
     4      pi* (taun4*(a29 + taun6*(a30 + taun4*a31)) +
     5      pi4*pi2*taun20*(taun9*a32 + taun20*taun10*(a33 +
     6         pi2*taun7*a34) +
     7      pi4*(a35 + taun15*(a36 + taun13*a37) +
     8      pi* (taun*a38 +
     9      pi* (taun19*(taun14*a39 +
     1          pi* a40) +
     2      pi2*(taun6*(a41 + taun14*(a42 +
     3      taun14*taun4*a43))) ))) ))) ))) ))) )
C
         gapt= e02 + taun*(e03 + taun*(e04 + taun3*e05)) +
     3      pi* (e06 + taun*(e07 + taun2*(e08 +
     4         taun3*(e09 + taun29*e10))) +
     5      pi* (e12 + taun2*(e13 + taun3*(e14 + taun29*e15)) +
     7      pi* (e16 + taun*(e17 + taun*e18) +
     8      pi* (taun6*e19 +
     9      pi* (taun2*(e20 + taun13*(e21 + taun19*e22)) +
     1      pi* (taun10*(e24 + taun14*e25) +
     2      pi* (taun6*(taun*e26 + taun29*e27) +
     3      pi* (taun6*taun6*e28 +
     4      pi* (taun3*(e29 + taun6*(e30 + taun4*e31)) +
     5      pi4*pi2*taun19*(taun9*e32 + taun20*taun10*(e33 +
     6         pi2*taun7*e34) +
     7      pi4*(e35 + taun15*(e36 + taun13*e37) +
     8      pi* (taun*e38 +
     9      pi* (taun19*(taun14*e39 + pi* e40) +
     2      pi2*(taun6*(e41 + taun14*(e42 +
     3      taun14*taun4*e43))) ))) ))) ))) ))) )
C
         GA0=1.D0
C
      ELSE
C
      TN=647.096D0
      DN=322.D0
C
      rho= d
      tau= tn/t 
      del= rho/DN
C
      tau2= tau*tau
      tau4= tau2*tau2
      tau12= tau4*tau4*tau4
C
      del2= del*del                                   
C
      phd= r01/del + tau2*(r09 + tau4*r10 + tau12*tau*(r11 +
     1      tau2*(r12 + del*tau4*tau*(r17 + tau4*(r18 + del*(r23 +
     2      del*(r27 + del*(r30 + del*(r33 + del2*(r35 +
     3      del*(r37 + del2*r40))) ))) ))) )) +
     4      del* (r13 + tau2*(r14 + tau4*(r15 + tau*r16)) +
     6      del* (r19 + tau2*(r20 + tau2*(r21 + tau12*r22)) +
     7      del* (r24 + tau2*(r25 + tau2*r26) +
     8      del* (tau*(r28 + tau2*r29) +
     1      del* (r31 + tau2*r32 +
     2      del* (tau2*r34 +
     3      del2*(tau2*r36 +
     4      del* (r38 + tau*r39))) ))) ))
C
      phdd= - s01/del2 + s13 + tau2*(s14 + tau4*(s15 + tau*s16 +
     1       tau12*tau4*(s17 + tau4*(s18 + del*(s23 + del*(s27 +
     2       del*(s30 + del*(s33 + del2*(s35 + del*(s37 +
     3       del2*s40))) ))) ))) ) +
     6       del* (s19 + tau2*(s20 + tau2*(s21 + tau12*s22)) +
     7       del* (s24 + tau2*(s25 + tau2*s26) +
     8       del* (tau*(s28 + tau2*s29) +
     1       del* (s31 + tau2*s32 +
     2       del* (tau2*s34 +
     3       del2*(tau2*s36 +
     4       del* (s38 + tau*s39))) ))) )
C
      phdt= tau*(t09 + tau4*t10 + tau12*tau*(t11 +
     1       tau2*(t12 + del*tau4*tau*(t17 + tau4*(t18 + del*(t23 +
     2       del*(t27 + del*(t30 + del*(t33 + del2*(t35 +
     3       del*(t37 + del2*t40))) ))) ))) )) +
     4       del* (tau*(t14 + tau4*(t15 + tau*t16)) +
     6       del* (tau*(t20 + tau2*(t21 + tau12*t22)) +
     7       del* (tau*(t25 + tau2*t26) +
     8       del* (t28 + tau2*t29 +
     1       del* (tau*t32 +
     2       del* (tau*t34 +
     3       del2*(tau*t36 +
     4       del*t39))) ))) )
C     
         DPDDEL=R*TN*DN/TAU*(2.D0*DEL*PHD+DEL2*PHDD)*1.D-6
         DPDD=DNR/DN/PNR*DPDDEL
C      
         DPDTAU=D*R*T/TAU*(DEL*TAU*PHDT-DEL*PHD)*1.D-6
         DPDT=-TN*TNR/PNR/T/T*DPDTAU
C
         GOTO 100
C
      ENDIF
C
      DVDPI=R*T/P/P*PN*(PI*PI*GAPP-GA0)*1.D-6
      DPDD=DNR/PNR/(-DVDPI*D*D/PN)
C
      DVDTAU=R*TN/PI/PN/TAU2*(GAPT*PI*TAU-PI*GAP-GA0)*1.D-6
      DDDT=TNR*TN/T/T/DNR*D*D*DVDTAU
C
      DPDT=DPDD*(-DDDT)
C
  100 CONTINUE
C  
      CHI=DR/DPDD
C      
      END
C
C************************************************************************
      DOUBLEPRECISION FUNCTION DIELC(T,D)
C************************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      PARAMETER (A01= 0.978224486826D0,
     1           A02=-0.957771379375D0,
     2           A03= 0.237511794148D0,
     3           A04= 0.714692244396D0,
     4           A05=-0.298217036956D0,
     5           A06=-0.108863472196D0,
     6           A07= 0.949327488264D-1,
     7           A08=-0.980469816509D-2,
     8           A09= 0.165167634970D-4,
     9           A10= 0.937359795772D-4,
     &           A11=-0.123179218720D-9,
     1           A12= 0.196096504426D-2)
c       
c     Die Koeffizienten zur Berechnung der Größen A und B wurden, soweit möglich, direkt
c     ausgerechnet. Die Berechnung der molaren Dichte war aufgrund dessen nicht mehr notwendig.
c      PARAMETER (CNA = 6.0221367D23,
c     1           CMDM= 6.138D-30,
c     2           ALP= 1.636D-40,
c     3           CK  = 1.380658D-23,
c     4           RI = 3.14159265359D0)

C
c      EPS0=1.D0/(4.D-7*RI*299792458.D0**2.D0)
C
      TR=647.096D0/T
      TR2=TR*TR
      TR4=TR2*TR2
      TRSR=DSQRT(TR)
      TRSR2=DSQRT(TRSR)
      DR=D/322.D0
C
      G=1.D0+DR*(A01*TRSR2+A02*TR+A03*TR2*TRSR+DR*(A04*TR*TRSR+
     1  DR*(TR*TRSR*(A05+TR*A06)+DR*(A07*TR2+DR*(A08*TR2+DR*(A09*
     2  TR4*TR+DR*(A10*TRSR+DR*DR*DR*A11*TR4*TR4*TR2)))))))+
     3  A12*DR/(T/228.D0-1.D0)**1.2D0
c
      A=1.592062421480592D-2*D*G*TR
      B=2.058842842420929D-4*D
C
      DIELC=(1.D0+A+5.D0*B+DSQRT(9.D0+2.D0*A+18.D0*B+A*A+
     1          10.D0*A*B+9.D0*B*B))/(4.D0-4.D0*B)
C
      END
C
C************************************************************************
      DOUBLEPRECISION FUNCTION REFI(T,D,W)
C************************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      PARAMETER (a0 = 0.244257733D0,
     1           a1 = 9.74634476D-3,
     2           a2 =-3.73234996D-3,
     3           a3 = 2.68678472D-4,
     4           a4 = 1.58920570D-3,
     5           a5 = 2.45934259D-3,
     6           a6 = 0.900704920D0,
     7           a7 =-1.66626219D-2)
      PARAMETER (CLUV= 0.229202D0,
     1           CLIR= 5.432937D0)
C
      TR=T/273.15D0
      DR=D/1000.D0
      WR=W/0.589D0
C
      WR2=WR*WR
C
      RIQ=DR*(a0+a1*DR+a2*TR+a3*WR2*TR+a4/WR2+a5/(WR2-CLUV*CLUV)+
     1        a6/(WR2-CLIR*CLIR)+a7*DR*DR)
      REFI=DSQRT((RIQ*2.D0+1.D0)/(1.D0-RIQ))
C
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION DVDPPT(P,V,T,IREG)
C************************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      PARAMETER (c15= -0.471 843 210 732 67 D-03*2.D0,
     6           c16= -0.300 017 807 930 26 D-03*2.D0,
     7           c17=  0.476 613 939 069 87 D-04*2.D0,
     8           c18= -0.441 418 453 308 46 D-05*2.D0,
     9           c19= -0.726 949 962 975 94 D-15*2.D0,
     1           c20= -0.316 796 448 450 54 D-04*6.D0,
     2           c21= -0.282 707 979 853 12 D-05*6.D0,
     3           c22= -0.852 051 281 201 03 D-09*6.D0,
     4           c23= -0.224 252 819 080 00 D-05*12.D0,
     5           c24= -0.651 712 228 956 01 D-06*12.D0,
     6           c25= -0.143 417 299 379 24 D-12*12.D0,
     7           c26= -0.405 169 968 601 17 D-06*20.D0,
     8           c27= -0.127 343 017 416 41 D-08*56.D0,
     9           c28= -0.174 248 712 306 34 D-09*56.D0,
     1           c29= -0.687 621 312 955 31 D-18*420.D0,
     2           c30=  0.144 783 078 285 21 D-19*506.D0,
     3           c31=  0.263 357 816 627 95 D-22*812.D0,
     4           c32= -0.119 476 226 400 71 D-22*870.D0,
     5           c33=  0.182 280 945 814 04 D-23*930.D0,
     6           c34= -0.935 370 872 924 58 D-25*992.D0)
C
      PARAMETER (a06=-0.330 326 416 702 03 D-04*2.D0,
     6           a07=-0.189 489 875 163 15 D-03*2.D0,
     7           a08=-0.393 927 772 433 55 D-02*2.D0,
     8           a09=-0.437 972 956 505 73 D-01*2.D0,
     9           a10=-0.266 745 479 140 87 D-04*2.D0,
     1           a11= 0.204 817 376 923 09 D-07*6.D0,
     2           a12= 0.438 706 672 844 35 D-06*6.D0,
     3           a13=-0.322 776 772 385 70 D-04*6.D0,
     4           a14=-0.150 339 245 421 48 D-02*6.D0,
     5           a15=-0.406 682 535 626 49 D-01*6.D0,
     6           a16=-0.788 473 095 593 67 D-09*12.D0,
     7           a17= 0.127 907 178 522 85 D-07*12.D0,
     8           a18= 0.482 253 727 185 07 D-06*12.D0,
     9           a19= 0.229 220 763 376 61 D-05*20.D0,
     1           a20=-0.167 147 664 510 61 D-10*30.D0,
     2           a21=-0.211 714 723 213 55 D-02*30.D0,
     3           a22=-0.238 957 419 341 04 D02*30.D0,
     4           a23=-0.590 595 643 242 70 D-17*42.D0,
     5           a24=-0.126 218 088 991 01 D-05*42.D0,
     6           a25=-0.389 468 424 357 39 D-01*42.D0,
     7           a26= 0.112 562 113 604 59 D-10*56.D0,
     8           a27=-0.823 113 408 979 98 D01*56.D0,
     9           a28= 0.198 097 128 020 88 D-07*72.D0,
     1           a29= 0.104 069 652 101 74 D-18*90.D0,
     2           a30=-0.102 347 470 959 29 D-12*90.D0,
     3           a31=-0.100 181 793 795 11 D-08*90.D0,
     4           a32=-0.808 829 086 469 85 D-10*240.D0,
     5           a33= 0.106 930 318 794 09 D0*240.D0,
     6           a34=-0.336 622 505 741 71 D0*306.D0,
     7           a35= 0.891 858 453 554 21 D-24*380.D0,
     8           a36= 0.306 293 168 762 32 D-12*380.D0,
     9           a37=-0.420 024 676 982 08 D-05*380.D0,
     1           a38=-0.590 560 296 856 39 D-25*420.D0,
     2           a39= 0.378 269 476 134 57 D-05*462.D0,
     3           a40=-0.127 686 089 346 81 D-14*506.D0,
     4           a41= 0.730 876 105 950 61 D-28*552.D0,
     5           a42= 0.554 147 153 507 78 D-16*552.D0,
     6           a43=-0.943 697 072 412 10 D-06*552.D0)
C
      PARAMETER (r01=  0.106 580 700 285 13 D01,
     1           r09= -0.126 543 154 777 14 D01,
     9           r10= -0.115 244 078 066 81 D01,
     1           r11=  0.885 210 439 843 18 D0,
     2           r12= -0.642 077 651 816 07 D0,
     3           r13=  0.384 934 601 866 71 D0*2.D0,
     4           r14= -0.852 147 088 242 06 D0*2.D0,
     5           r15=  0.489 722 815 418 77 D01*2.D0,
     6           r16= -0.305 026 172 569 65 D01*2.D0,
     7           r17=  0.394 205 368 791 54 D-01*2.D0,
     8           r18=  0.125 584 084 243 08 D0*2.D0,
     9           r19= -0.279 993 296 987 10 D0*3.D0,
     1           r20=  0.138 997 995 694 60 D01*3.D0,
     2           r21= -0.201 899 150 235 70 D01*3.D0,
     3           r22= -0.821 476 371 739 63 D-02*3.D0,
     4           r23= -0.475 960 357 349 23 D0*3.D0,
     5           r24=  0.439 840 744 735 00 D-01*4.D0,
     6           r25= -0.444 764 354 287 39 D0*4.D0,
     7           r26=  0.905 720 707 197 33 D0*4.D0,
     8           r27=  0.705 224 500 879 67 D0*4.D0,
     9           r28=  0.107 705 126 263 32 D0*5.D0,
     1           r29= -0.329 136 232 589 54 D0*5.D0,
     2           r30= -0.508 710 620 411 58 D0*5.D0,
     3           r31= -0.221 754 008 730 96 D-01*6.D0,
     4           r32=  0.942 607 516 650 92 D-01*6.D0,
     5           r33=  0.164 362 784 479 61 D0*6.D0,
     6           r34= -0.135 033 722 413 48 D-01*7.D0,
     7           r35= -0.148 343 453 524 72 D-01*8.D0,
     8           r36=  0.579 229 536 280 84 D-03*9.D0,
     9           r37=  0.323 089 047 037 11 D-02*9.D0,
     1           r38=  0.809 648 029 962 15 D-04*10.D0,
     2           r39= -0.165 576 797 950 37 D-03*10.D0,
     3           r40= -0.449 238 990 618 15 D-04*11.D0)
      PARAMETER (s01=  0.106 580 700 285 13 D01,
     3           s13=  0.384 934 601 866 71 D0*2.D0,
     4           s14= -0.852 147 088 242 06 D0*2.D0,
     5           s15=  0.489 722 815 418 77 D01*2.D0,
     6           s16= -0.305 026 172 569 65 D01*2.D0,
     7           s17=  0.394 205 368 791 54 D-01*2.D0,
     8           s18=  0.125 584 084 243 08 D0*2.D0,
     9           s19= -0.279 993 296 987 10 D0*6.D0,
     1           s20=  0.138 997 995 694 60 D01*6.D0,
     2           s21= -0.201 899 150 235 70 D01*6.D0,
     3           s22= -0.821 476 371 739 63 D-02*6.D0,
     4           s23= -0.475 960 357 349 23 D0*6.D0,
     5           s24=  0.439 840 744 735 00 D-01*12.D0,
     6           s25= -0.444 764 354 287 39 D0*12.D0,
     7           s26=  0.905 720 707 197 33 D0*12.D0,
     8           s27=  0.705 224 500 879 67 D0*12.D0,
     9           s28=  0.107 705 126 263 32 D0*20.D0,
     1           s29= -0.329 136 232 589 54 D0*20.D0,
     2           s30= -0.508 710 620 411 58 D0*20.D0,
     3           s31= -0.221 754 008 730 96 D-01*30.D0,
     4           s32=  0.942 607 516 650 92 D-01*30.D0,
     5           s33=  0.164 362 784 479 61 D0*30.D0,
     6           s34= -0.135 033 722 413 48 D-01*42.D0,
     7           s35= -0.148 343 453 524 72 D-01*56.D0,
     8           s36=  0.579 229 536 280 84 D-03*72.D0,
     9           s37=  0.323 089 047 037 11 D-02*72.D0,
     1           s38=  0.809 648 029 962 15 D-04*90.D0,
     2           s39= -0.165 576 797 950 37 D-03*90.D0,
     3           s40= -0.449 238 990 618 15 D-04*110.D0)
C
      PARAMETER (H10= -.36668082266957D-05,
     &           H11= .17887679267013D-06, R1=461.526D-6)
C
      PARAMETER (R=0.461526D3)
C
      IF (IREG .EQ. 1) THEN
C
         TN= 1386.D0
         PN= 16.53D0
C
         tau= tn/t
         TAU2=TAU*TAU
         taun= tau - 1.222d0
         taun2= taun*taun
         taun3= taun2*taun
         taun6= taun3*taun3
         taun12= taun6*taun6
         tauninv= 1.d0/taun
         tauninv2= tauninv*tauninv
         tauninv3= tauninv2*tauninv
         tauninv7= tauninv3*tauninv3*tauninv
C
         PI=P/PN
         pin= 7.1d0 - p/pn
         pin2= pin*pin
         pin3= pin2*pin
         pin6= pin3*pin3
C         
         gapp= tauninv3*(c15 + taun3*(c16 + taun*(c17 +
     5       taun2*(c18 + taun12*taun2*c19))) +
     6       pin*  tauninv*(c20 + taun3*taun*(c21 + taun6*c22) +
     7       pin*  tauninv*(c23 + taun3*(c24 + taun12*c25) +
     8       pin*  tauninv3*(c26 +
     9       pin3* (tauninv3*(c27 + taun3*taun2*c28) +
     1       pin6*pin6*pin*tauninv7*tauninv7*tauninv7*(c29 +
     2       pin2* tauninv2*(c30 +
     3       pin6* tauninv7*(c31 +
     4       pin*  tauninv*(c32 +
     5       pin*  tauninv*(c33 +
     6       pin*  tauninv*c34))) ))) ))) )
C
         GA0=0.D0
C
      ELSE IF (IREG .EQ. 2) THEN
C
         TN= 540.D0
         PN= 1.D0
C
         pi= p/pn
         tau=tn/t
         TAU2=TAU*TAU
C
         taun= tau - 0.5d0
         taun2= taun*taun
         taun3= taun2*taun
         taun4= taun3*taun
         taun6= taun4*taun2
         taun7= taun6*taun
         taun9= taun7*taun2
         taun10= taun9*taun
         taun13= taun10*taun3
         taun14= taun13*taun
         taun15= taun14*taun
         taun19= taun15*taun4
         taun20= taun19*taun
         taun29= taun20*taun9
C
         pi2= pi*pi
         pi4= pi2*pi2
C
         gapp= (taun*(a06 + taun*(a07 + taun2*(a08 +
     4      taun3*(a09 + taun29*a10)))) +
     5      pi* (a11 + taun*(a12 + taun2*(a13 + taun3*(a14 +
     6         taun29*a15))) +
     7      pi* (taun*(a16 + taun*(a17 + taun*a18)) +
     8      pi* (taun7*a19 +
     9      pi* (taun3*(a20 + taun13*(a21 + taun19*a22)) +
     1      pi* (a23 + taun10*taun*(a24 + taun14*a25) +
     2      pi* (taun7*(taun*a26 + taun29*a27) +
     3      pi* (taun13*a28 +
     4      pi* (taun4*(a29 + taun6*(a30 + taun4*a31)) +
     5      pi4*pi2*taun20*(taun9*a32 + taun20*taun10*(a33 +
     6         pi2*taun7*a34) +
     7      pi4*(a35 + taun15*(a36 + taun13*a37) +
     8      pi* (taun*a38 +
     9      pi* (taun19*(taun14*a39 +
     1          pi* a40) +
     2      pi2*(taun6*(a41 + taun14*(a42 +
     3      taun14*taun4*a43))) ))) ))) ))) ))) )
C
         GA0=1.D0
C
      ELSEIF (IREG .EQ. 3) THEN
C
      TN=647.096D0
      DN=322.D0
C
      rho= 1.d0/v
      tau= tn/t 
      del= rho/DN
C
      tau2= tau*tau
      tau4= tau2*tau2
      tau12= tau4*tau4*tau4
C
      del2= del*del                                   
C
      phd= r01/del + tau2*(r09 + tau4*r10 + tau12*tau*(r11 +
     1      tau2*(r12 + del*tau4*tau*(r17 + tau4*(r18 + del*(r23 +
     2      del*(r27 + del*(r30 + del*(r33 + del2*(r35 +
     3      del*(r37 + del2*r40))) ))) ))) )) +
     4      del* (r13 + tau2*(r14 + tau4*(r15 + tau*r16)) +
     6      del* (r19 + tau2*(r20 + tau2*(r21 + tau12*r22)) +
     7      del* (r24 + tau2*(r25 + tau2*r26) +
     8      del* (tau*(r28 + tau2*r29) +
     1      del* (r31 + tau2*r32 +
     2      del* (tau2*r34 +
     3      del2*(tau2*r36 +
     4      del* (r38 + tau*r39))) ))) ))
C
      phdd= - s01/del2 + s13 + tau2*(s14 + tau4*(s15 + tau*s16 +
     1       tau12*tau4*(s17 + tau4*(s18 + del*(s23 + del*(s27 +
     2       del*(s30 + del*(s33 + del2*(s35 + del*(s37 +
     3       del2*s40))) ))) ))) ) +
     6       del* (s19 + tau2*(s20 + tau2*(s21 + tau12*s22)) +
     7       del* (s24 + tau2*(s25 + tau2*s26) +
     8       del* (tau*(s28 + tau2*s29) +
     1       del* (s31 + tau2*s32 +
     2       del* (tau2*s34 +
     3       del2*(tau2*s36 +
     4       del* (s38 + tau*s39))) ))) )
C
         DPDDEL=R*TN*DN/TAU*(2.D0*DEL*PHD+DEL2*PHDD)
         DVDPPT=1.D0/DN/DEL2/(-DPDDEL)*1.D6
C
c dvdp in m^3/(kg*MPa)
C
         GOTO 100
C
      ELSEIF (IREG .EQ. 5) THEN
C
      TAU=T*1.D-3
      TAUM2=1.D0/(TAU*TAU)
C
      DVDPPT=(H10*TAUM2*TAUM2*TAUM2+2.D0*H11*P)*TAUM2-R1*T/P/P
C
c dvdp in m^3/(kg*MPa)
C
         GOTO 100
C
      ENDIF
C
      DVDPI=R*T/P/P*PN*(PI*PI*GAPP-GA0)
      DVDPPT=DVDPI/PN*1.D-6
C
c dvdp in m^3/(kg*MPa)
C
  100 CONTINUE
C  
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION DVDTPT(P,V,T,IREG)
C************************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      PARAMETER (b09=  0.283 190 801 238 04 D-03,
     9           b10= -0.607 063 015 658 74 D-03,
     1           b11= -0.189 900 682 184 19 D-01,
     2           b12= -0.325 297 487 705 05 D-01,
     3           b13= -0.218 417 171 754 14 D-01,
     4           b14= -0.528 383 579 699 30 D-04,
     5           b15= -0.471 843 210 732 67 D-03*2.D0,
     6           b16= -0.300 017 807 930 26 D-03*2.D0,
     7           b17=  0.476 613 939 069 87 D-04*2.D0,
     8           b18= -0.441 418 453 308 46 D-05*2.D0,
     9           b19= -0.726 949 962 975 94 D-15*2.D0,
     1           b20= -0.316 796 448 450 54 D-04*3.D0,
     2           b21= -0.282 707 979 853 12 D-05*3.D0,
     3           b22= -0.852 051 281 201 03 D-09*3.D0,
     4           b23= -0.224 252 819 080 00 D-05*4.D0,
     5           b24= -0.651 712 228 956 01 D-06*4.D0,
     6           b25= -0.143 417 299 379 24 D-12*4.D0,
     7           b26= -0.405 169 968 601 17 D-06*5.D0,
     8           b27= -0.127 343 017 416 41 D-08*8.D0,
     9           b28= -0.174 248 712 306 34 D-09*8.D0,
     1           b29= -0.687 621 312 955 31 D-18*21.D0,
     2           b30=  0.144 783 078 285 21 D-19*23.D0,
     3           b31=  0.263 357 816 627 95 D-22*29.D0,
     4           b32= -0.119 476 226 400 71 D-22*30.D0,
     5           b33=  0.182 280 945 814 04 D-23*31.D0,
     6           b34= -0.935 370 872 924 58 D-25*32.D0)
C
      PARAMETER (f09= -0.283 190 801 238 04 D-03*9.D0,
     9           f10=  0.607 063 015 658 74 D-03*7.D0,
     1           f11=  0.189 900 682 184 19 D-01*1.D0,
     3           f13= -0.218 417 171 754 14 D-01,
     4           f14= -0.528 383 579 699 30 D-04*3.D0,
     5           f15=  0.471 843 210 732 67 D-03*6.D0,
     7           f17=  0.476 613 939 069 87 D-04*2.D0,
     8           f18= -0.441 418 453 308 46 D-05*6.D0,
     9           f19= -0.726 949 962 975 94 D-15*34.D0,
     1           f20=  0.316 796 448 450 54 D-04*12.D0,
     3           f22= -0.852 051 281 201 03 D-09*18.D0,
     4           f23=  0.224 252 819 080 00 D-05*20.D0,
     5           f24=  0.651 712 228 956 01 D-06*8.D0,
     6           f25= -0.143 417 299 379 24 D-12*40.D0,
     7           f26=  0.405 169 968 601 17 D-06*40.D0,
     8           f27=  0.127 343 017 416 41 D-08*88.D0,
     9           f28=  0.174 248 712 306 34 D-09*48.D0,
     1           f29=  0.687 621 312 955 31 D-18*609.D0,
     2           f30= -0.144 783 078 285 21 D-19*713.D0,
     3           f31= -0.263 357 816 627 95 D-22*1102.D0,
     4           f32=  0.119 476 226 400 71 D-22*1170.D0,
     5           f33= -0.182 280 945 814 04 D-23*1240.D0,
     6           f34=  0.935 370 872 924 58 D-25*1312.D0)
C
      PARAMETER (e02=-0.178 348 622 923 58 D-01,
     2           e03=-0.459 960 136 963 65 D-01*2.D0,
     3           e04=-0.575 812 590 834 32 D-01*3.D0,
     4           e05=-0.503 252 787 279 30 D-01*6.D0,
     5           e06=-0.330 326 416 702 03 D-04*2.D0,
     6           e07=-0.189 489 875 163 15 D-03*4.D0,
     7           e08=-0.393 927 772 433 55 D-02*8.D0,
     8           e09=-0.437 972 956 505 73 D-01*14.D0,
     9           e10=-0.266 745 479 140 87 D-04*72.D0,
     2           e12= 0.438 706 672 844 35 D-06*3.D0,
     3           e13=-0.322 776 772 385 70 D-04*9.D0,
     4           e14=-0.150 339 245 421 48 D-02*18.D0,
     5           e15=-0.406 682 535 626 49 D-01*105.D0,
     6           e16=-0.788 473 095 593 67 D-09*4.D0,
     7           e17= 0.127 907 178 522 85 D-07*8.D0,
     8           e18= 0.482 253 727 185 07 D-06*12.D0,
     9           e19= 0.229 220 763 376 61 D-05*35.D0,
     1           e20=-0.167 147 664 510 61 D-10*18.D0,
     2           e21=-0.211 714 723 213 55 D-02*96.D0,
     3           e22=-0.238 957 419 341 04 D02*210.D0,
     5           e24=-0.126 218 088 991 01 D-05*77.D0,
     6           e25=-0.389 468 424 357 39 D-01*175.D0,
     7           e26= 0.112 562 113 604 59 D-10*64.D0,
     8           e27=-0.823 113 408 979 98 D01*288.D0,
     9           e28= 0.198 097 128 020 88 D-07*117.D0,
     1           e29= 0.104 069 652 101 74 D-18*40.D0,
     2           e30=-0.102 347 470 959 29 D-12*100.D0,
     3           e31=-0.100 181 793 795 11 D-08*140.D0,
     4           e32=-0.808 829 086 469 85 D-10*464.D0,
     5           e33= 0.106 930 318 794 09 D0*800.D0,
     6           e34=-0.336 622 505 741 71 D0*1026.D0,
     7           e35= 0.891 858 453 554 21 D-24*400.D0,
     8           e36= 0.306 293 168 762 32 D-12*700.D0,
     9           e37=-0.420 024 676 982 08 D-05*960.D0,
     1           e38=-0.590 560 296 856 39 D-25*441.D0,
     2           e39= 0.378 269 476 134 57 D-05*1166.D0,
     3           e40=-0.127 686 089 346 81 D-14*897.D0,
     4           e41= 0.730 876 105 950 61 D-28*624.D0,
     5           e42= 0.554 147 153 507 78 D-16*960.D0,
     6           e43=-0.943 697 072 412 10 D-06*1392.D0)
C
      PARAMETER (d01=-0.177 317 424 732 13 D-02,
     1           d02=-0.178 348 622 923 58 D-01,
     2           d03=-0.459 960 136 963 65 D-01,
     3           d04=-0.575 812 590 834 32 D-01,
     4           d05=-0.503 252 787 279 30 D-01,
     5           d06=-0.330 326 416 702 03 D-04*2.D0,
     6           d07=-0.189 489 875 163 15 D-03*2.D0,
     7           d08=-0.393 927 772 433 55 D-02*2.D0,
     8           d09=-0.437 972 956 505 73 D-01*2.D0,
     9           d10=-0.266 745 479 140 87 D-04*2.D0,
     1           d11= 0.204 817 376 923 09 D-07*3.D0,
     2           d12= 0.438 706 672 844 35 D-06*3.D0,
     3           d13=-0.322 776 772 385 70 D-04*3.D0,
     4           d14=-0.150 339 245 421 48 D-02*3.D0,
     5           d15=-0.406 682 535 626 49 D-01*3.D0,
     6           d16=-0.788 473 095 593 67 D-09*4.D0,
     7           d17= 0.127 907 178 522 85 D-07*4.D0,
     8           d18= 0.482 253 727 185 07 D-06*4.D0,
     9           d19= 0.229 220 763 376 61 D-05*5.D0,
     1           d20=-0.167 147 664 510 61 D-10*6.D0,
     2           d21=-0.211 714 723 213 55 D-02*6.D0,
     3           d22=-0.238 957 419 341 04 D02*6.D0,
     4           d23=-0.590 595 643 242 70 D-17*7.D0,
     5           d24=-0.126 218 088 991 01 D-05*7.D0,
     6           d25=-0.389 468 424 357 39 D-01*7.D0,
     7           d26= 0.112 562 113 604 59 D-10*8.D0,
     8           d27=-0.823 113 408 979 98 D01*8.D0,
     9           d28= 0.198 097 128 020 88 D-07*9.D0,
     1           d29= 0.104 069 652 101 74 D-18*10.D0,
     2           d30=-0.102 347 470 959 29 D-12*10.D0,
     3           d31=-0.100 181 793 795 11 D-08*10.D0,
     4           d32=-0.808 829 086 469 85 D-10*16.D0,
     5           d33= 0.106 930 318 794 09 D0*16.D0,
     6           d34=-0.336 622 505 741 71 D0*18.D0,
     7           d35= 0.891 858 453 554 21 D-24*20.D0,
     8           d36= 0.306 293 168 762 32 D-12*20.D0,
     9           d37=-0.420 024 676 982 08 D-05*20.D0,
     1           d38=-0.590 560 296 856 39 D-25*21.D0,
     2           d39= 0.378 269 476 134 57 D-05*22.D0,
     3           d40=-0.127 686 089 346 81 D-14*23.D0,
     4           d41= 0.730 876 105 950 61 D-28*24.D0,
     5           d42= 0.554 147 153 507 78 D-16*24.D0,
     6           d43=-0.943 697 072 412 10 D-06*24.D0)
C
      PARAMETER (r01=  0.106 580 700 285 13 D01,
     1           r09= -0.126 543 154 777 14 D01,
     9           r10= -0.115 244 078 066 81 D01,
     1           r11=  0.885 210 439 843 18 D0,
     2           r12= -0.642 077 651 816 07 D0,
     3           r13=  0.384 934 601 866 71 D0*2.D0,
     4           r14= -0.852 147 088 242 06 D0*2.D0,
     5           r15=  0.489 722 815 418 77 D01*2.D0,
     6           r16= -0.305 026 172 569 65 D01*2.D0,
     7           r17=  0.394 205 368 791 54 D-01*2.D0,
     8           r18=  0.125 584 084 243 08 D0*2.D0,
     9           r19= -0.279 993 296 987 10 D0*3.D0,
     1           r20=  0.138 997 995 694 60 D01*3.D0,
     2           r21= -0.201 899 150 235 70 D01*3.D0,
     3           r22= -0.821 476 371 739 63 D-02*3.D0,
     4           r23= -0.475 960 357 349 23 D0*3.D0,
     5           r24=  0.439 840 744 735 00 D-01*4.D0,
     6           r25= -0.444 764 354 287 39 D0*4.D0,
     7           r26=  0.905 720 707 197 33 D0*4.D0,
     8           r27=  0.705 224 500 879 67 D0*4.D0,
     9           r28=  0.107 705 126 263 32 D0*5.D0,
     1           r29= -0.329 136 232 589 54 D0*5.D0,
     2           r30= -0.508 710 620 411 58 D0*5.D0,
     3           r31= -0.221 754 008 730 96 D-01*6.D0,
     4           r32=  0.942 607 516 650 92 D-01*6.D0,
     5           r33=  0.164 362 784 479 61 D0*6.D0,
     6           r34= -0.135 033 722 413 48 D-01*7.D0,
     7           r35= -0.148 343 453 524 72 D-01*8.D0,
     8           r36=  0.579 229 536 280 84 D-03*9.D0,
     9           r37=  0.323 089 047 037 11 D-02*9.D0,
     1           r38=  0.809 648 029 962 15 D-04*10.D0,
     2           r39= -0.165 576 797 950 37 D-03*10.D0,
     3           r40= -0.449 238 990 618 15 D-04*11.D0)
      PARAMETER (s01=  0.106 580 700 285 13 D01,
     3           s13=  0.384 934 601 866 71 D0*2.D0,
     4           s14= -0.852 147 088 242 06 D0*2.D0,
     5           s15=  0.489 722 815 418 77 D01*2.D0,
     6           s16= -0.305 026 172 569 65 D01*2.D0,
     7           s17=  0.394 205 368 791 54 D-01*2.D0,
     8           s18=  0.125 584 084 243 08 D0*2.D0,
     9           s19= -0.279 993 296 987 10 D0*6.D0,
     1           s20=  0.138 997 995 694 60 D01*6.D0,
     2           s21= -0.201 899 150 235 70 D01*6.D0,
     3           s22= -0.821 476 371 739 63 D-02*6.D0,
     4           s23= -0.475 960 357 349 23 D0*6.D0,
     5           s24=  0.439 840 744 735 00 D-01*12.D0,
     6           s25= -0.444 764 354 287 39 D0*12.D0,
     7           s26=  0.905 720 707 197 33 D0*12.D0,
     8           s27=  0.705 224 500 879 67 D0*12.D0,
     9           s28=  0.107 705 126 263 32 D0*20.D0,
     1           s29= -0.329 136 232 589 54 D0*20.D0,
     2           s30= -0.508 710 620 411 58 D0*20.D0,
     3           s31= -0.221 754 008 730 96 D-01*30.D0,
     4           s32=  0.942 607 516 650 92 D-01*30.D0,
     5           s33=  0.164 362 784 479 61 D0*30.D0,
     6           s34= -0.135 033 722 413 48 D-01*42.D0,
     7           s35= -0.148 343 453 524 72 D-01*56.D0,
     8           s36=  0.579 229 536 280 84 D-03*72.D0,
     9           s37=  0.323 089 047 037 11 D-02*72.D0,
     1           s38=  0.809 648 029 962 15 D-04*90.D0,
     2           s39= -0.165 576 797 950 37 D-03*90.D0,
     3           s40= -0.449 238 990 618 15 D-04*110.D0)
      PARAMETER (t09= -0.126 543 154 777 14 D01*2.D0,
     9           t10= -0.115 244 078 066 81 D01*6.D0,
     1           t11=  0.885 210 439 843 18 D0*15.D0,
     2           t12= -0.642 077 651 816 07 D0*17.D0,
     4           t14= -0.852 147 088 242 06 D0*4.D0,
     5           t15=  0.489 722 815 418 77 D01*12.D0,
     6           t16= -0.305 026 172 569 65 D01*14.D0,
     7           t17=  0.394 205 368 791 54 D-01*44.D0,
     8           t18=  0.125 584 084 243 08 D0*52.D0,
     1           t20=  0.138 997 995 694 60 D01*6.D0,
     2           t21= -0.201 899 150 235 70 D01*12.D0,
     3           t22= -0.821 476 371 739 63 D-02*48.D0,
     4           t23= -0.475 960 357 349 23 D0*78.D0,
     6           t25= -0.444 764 354 287 39 D0*8.D0,
     7           t26=  0.905 720 707 197 33 D0*16.D0,
     8           t27=  0.705 224 500 879 67 D0*104.D0,
     9           t28=  0.107 705 126 263 32 D0*5.D0,
     1           t29= -0.329 136 232 589 54 D0*15.D0,
     2           t30= -0.508 710 620 411 58 D0*130.D0,
     4           t32=  0.942 607 516 650 92 D-01*12.D0,
     5           t33=  0.164 362 784 479 61 D0*156.D0,
     6           t34= -0.135 033 722 413 48 D-01*14.D0,
     7           t35= -0.148 343 453 524 72 D-01*208.D0,
     8           t36=  0.579 229 536 280 84 D-03*18.D0,
     9           t37=  0.323 089 047 037 11 D-02*234.D0,
     2           t39= -0.165 576 797 950 37 D-03*10.D0,
     3           t40= -0.449 238 990 618 15 D-04*286.D0)
C
      PARAMETER (H7= -.21203806358652D-02,
     &           H9= -.57982358693700D-04, H10= -.36668082266957D-05,
     &           H11= .17887679267013D-06, R1=461.526D-6)
C
      PARAMETER (R=0.461526D3)
C
      IF (IREG .EQ. 1) THEN
C
         TN= 1386.D0
         PN= 16.53D0
C
         tau= tn/t
         TAU2=TAU*TAU
         taun= tau - 1.222d0
         taun2= taun*taun
         taun3= taun2*taun
         taun6= taun3*taun3
         taun12= taun6*taun6
         tauninv= 1.d0/taun
         tauninv2= tauninv*tauninv
         tauninv3= tauninv2*tauninv
         tauninv7= tauninv3*tauninv3*tauninv
C
         PI=P/PN
         pin= 7.1d0 - p/pn
         pin2= pin*pin
         pin3= pin2*pin
         pin6= pin3*pin3
C
         gap= -(tauninv7*(tauninv2*b09 + b10 + taun6*(b11 +
     1      taun*(b12 + taun*(b13 + taun2*b14)))) +
     2      pin*  tauninv3*(b15 + taun3*(b16 + taun*(b17 +
     3            taun2*(b18 + taun12*taun2*b19))) +
     4      pin*  tauninv*(b20 + taun3*taun*(b21 + taun6*b22) +
     5      pin*  tauninv*(b23 + taun3*(b24 + taun12*b25) +
     6      pin*  tauninv3*(b26 +
     7      pin3* (tauninv3*(b27 + taun3*taun2*b28) +
     8      pin6*pin6*pin*tauninv7*tauninv7*tauninv7*(b29 +
     9      pin2* tauninv2*(b30 +
     1      pin6* tauninv7*(b31 +
     2      pin*  tauninv*(b32 +
     3      pin*  tauninv*(b33 +
     4      pin*  tauninv*b34))) ))) ))) ))
C
         gapt= -(tauninv*(tauninv7*(tauninv2*f09 + f10 + taun6*(f11 +
     3       taun2*(f13 + taun2*f14))) +
     4       pin*  tauninv3*(f15 + taun3*taun*(f17 +
     5             taun2*(f18 + taun12*taun2*f19)) +
     6       pin*  tauninv*(f20 + taun6*taun3*taun*f22 +
     7       pin*  tauninv*(f23 + taun3*(f24 + taun12*f25) +
     8       pin*  tauninv3*(f26 +
     9       pin3*(tauninv3*(f27 + taun3*taun2*f28) +
     1       pin6*pin6*pin*tauninv7*tauninv7*tauninv7*(f29 +
     2       pin2* tauninv2*(f30 +
     3       pin6* tauninv7*(f31 +
     4       pin*  tauninv*(f32 +
     5       pin*  tauninv*(f33 +
     6       pin*  tauninv*f34))) ))) ))) )))
C     
         GA0=0.D0
C
      ELSE IF (IREG .EQ. 2) THEN
C
         TN= 540.D0
         PN= 1.D0
C
         pi= p/pn
         tau=tn/t
         TAU2=TAU*TAU
C
         taun= tau - 0.5d0
         taun2= taun*taun
         taun3= taun2*taun
         taun4= taun3*taun
         taun6= taun4*taun2
         taun7= taun6*taun
         taun9= taun7*taun2
         taun10= taun9*taun
         taun13= taun10*taun3
         taun14= taun13*taun
         taun15= taun14*taun
         taun19= taun15*taun4
         taun20= taun19*taun
         taun29= taun20*taun9
C
         pi2= pi*pi
         pi4= pi2*pi2
C
         gap= d01 + taun*(d02 + taun*(d03 + taun*(d04 + taun3*d05))) +
     3      pi* (taun*(d06 + taun*(d07 + taun2*(d08 +
     4         taun3*(d09 + taun29*d10)))) +
     5      pi* (d11 + taun*(d12 + taun2*(d13 + taun3*(d14 +
     6         taun29*d15))) +
     7      pi* (taun*(d16 + taun*(d17 + taun*d18)) +
     8      pi* (taun7*d19 +
     9      pi* (taun3*(d20 + taun13*(d21 + taun19*d22)) +
     1      pi* (d23 + taun10*taun*(d24 + taun14*d25) +
     2      pi* (taun7*(taun*d26 + taun29*d27) +
     3      pi* (taun13*d28 +
     4      pi* (taun4*(d29 + taun6*(d30 + taun4*d31)) +
     5      pi4*pi2*taun20*(taun9*d32 + taun20*taun10*(d33 +
     6         pi2*taun7*d34) +
     7      pi4*(d35 + taun15*(d36 + taun13*d37) +
     8      pi* (taun*d38 +
     9      pi* (taun19*(taun14*d39 +
     1          pi* d40) +
     2      pi2*(taun6*(d41 + taun14*(d42 +
     3      taun14*taun4*d43))) ))) ))) ))) ))) )
C
         gapt= e02 + taun*(e03 + taun*(e04 + taun3*e05)) +
     3      pi* (e06 + taun*(e07 + taun2*(e08 +
     4         taun3*(e09 + taun29*e10))) +
     5      pi* (e12 + taun2*(e13 + taun3*(e14 + taun29*e15)) +
     7      pi* (e16 + taun*(e17 + taun*e18) +
     8      pi* (taun6*e19 +
     9      pi* (taun2*(e20 + taun13*(e21 + taun19*e22)) +
     1      pi* (taun10*(e24 + taun14*e25) +
     2      pi* (taun6*(taun*e26 + taun29*e27) +
     3      pi* (taun6*taun6*e28 +
     4      pi* (taun3*(e29 + taun6*(e30 + taun4*e31)) +
     5      pi4*pi2*taun19*(taun9*e32 + taun20*taun10*(e33 +
     6         pi2*taun7*e34) +
     7      pi4*(e35 + taun15*(e36 + taun13*e37) +
     8      pi* (taun*e38 +
     9      pi* (taun19*(taun14*e39 + pi* e40) +
     2      pi2*(taun6*(e41 + taun14*(e42 +
     3      taun14*taun4*e43))) ))) ))) ))) ))) )
C
         GA0=1.D0
C
      ELSEIF (IREG .EQ. 3) THEN
C
      TN=647.096D0
      DN=322.D0
C
      d=1.d0/v
      rho= d
      tau= tn/t 
      del= rho/DN
C
      tau2= tau*tau
      tau4= tau2*tau2
      tau12= tau4*tau4*tau4
C
      del2= del*del                                   
C
      phd= r01/del + tau2*(r09 + tau4*r10 + tau12*tau*(r11 +
     1      tau2*(r12 + del*tau4*tau*(r17 + tau4*(r18 + del*(r23 +
     2      del*(r27 + del*(r30 + del*(r33 + del2*(r35 +
     3      del*(r37 + del2*r40))) ))) ))) )) +
     4      del* (r13 + tau2*(r14 + tau4*(r15 + tau*r16)) +
     6      del* (r19 + tau2*(r20 + tau2*(r21 + tau12*r22)) +
     7      del* (r24 + tau2*(r25 + tau2*r26) +
     8      del* (tau*(r28 + tau2*r29) +
     1      del* (r31 + tau2*r32 +
     2      del* (tau2*r34 +
     3      del2*(tau2*r36 +
     4      del* (r38 + tau*r39))) ))) ))
C
      phdd= - s01/del2 + s13 + tau2*(s14 + tau4*(s15 + tau*s16 +
     1       tau12*tau4*(s17 + tau4*(s18 + del*(s23 + del*(s27 +
     2       del*(s30 + del*(s33 + del2*(s35 + del*(s37 +
     3       del2*s40))) ))) ))) ) +
     6       del* (s19 + tau2*(s20 + tau2*(s21 + tau12*s22)) +
     7       del* (s24 + tau2*(s25 + tau2*s26) +
     8       del* (tau*(s28 + tau2*s29) +
     1       del* (s31 + tau2*s32 +
     2       del* (tau2*s34 +
     3       del2*(tau2*s36 +
     4       del* (s38 + tau*s39))) ))) )
C
      phdt= tau*(t09 + tau4*t10 + tau12*tau*(t11 +
     1       tau2*(t12 + del*tau4*tau*(t17 + tau4*(t18 + del*(t23 +
     2       del*(t27 + del*(t30 + del*(t33 + del2*(t35 +
     3       del*(t37 + del2*t40))) ))) ))) )) +
     4       del* (tau*(t14 + tau4*(t15 + tau*t16)) +
     6       del* (tau*(t20 + tau2*(t21 + tau12*t22)) +
     7       del* (tau*(t25 + tau2*t26) +
     8       del* (t28 + tau2*t29 +
     1       del* (tau*t32 +
     2       del* (tau*t34 +
     3       del2*(tau*t36 +
     4       del*t39))) ))) )
C     
         DPDDEL=R*TN*DN/TAU*(2.D0*DEL*PHD+DEL2*PHDD)
         DPDV=DPDDEL*(-DEL)*D
C      
c dpdv in Pa*kg/m^3
C
         DPDTAU=D*R*T/TAU*(DEL*TAU*PHDT-DEL*PHD)
         DPDT=DPDTAU*(-TAU)/T
C
c dpdt in Pa/K
C
         DVDTPT=(-DPDT)/DPDV
C
c dvdt in m^3/(kg*K)
C
         GOTO 100
C
      ELSEIF (IREG .EQ. 5) THEN
C
      TAU=T*1.D-3
      TAUM=1.D0/TAU
      TAUM2=TAUM*TAUM
C
      DVDTPT=-((H10*TAUM2*TAUM2*TAUM2*4.D0+H11*P)*P+H7)*TAUM2*TAUM*2.D-3
     &      +R1/P+H9*1.D-3
C
c dvdt in m^3/(kg*K)
C
         GOTO 100
C
      ENDIF
C
      DVDTAU=R*TN/PI/(PN*1.D6)/TAU2*(GAPT*PI*TAU-PI*GAP-GA0)
      DVDTPT=DVDTAU*(-TAU)/T
C
c dvdt in m^3/(kg*K)
C
  100 CONTINUE
C  
      END
C
C************************************************************************
      DOUBLEPRECISION FUNCTION SURFTT(T)
C************************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      TAU = (1.D0 - T/647.096D0)
      SURFTT = 235.8D0 *TAU ** 1.256D0 * (1.D0 - 0.625D+0 * TAU)
C
      END
c
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is pressure in MPa.
C  Second argument is temperature in K.
C  Returns specific volume in m^3/kg
C
      FUNCTION vpt1n(P,T)
C
      IMPLICIT double precision (A-H,O-Z)
C
      parameter(a8=-0.79068310787773D-08,
     9          a9= 0.16949507886565D-07,
     &         a10= 0.53021235478366D-06,
     &         a11= 0.90824711621636D-06,
     2         a12= 0.60983184277677D-06,
     3         a13= 0.14752738052287D-08,
     4         a14= 0.26348204437581D-07,
     &         a15= 0.16753299313106D-07,
     6         a16=-0.26614606756583D-08,
     7         a17= 0.24649255061299D-09,
     8         a18= 0.40593624756496D-19,
     9         a19= 0.26535353478691D-08,
     &         a20= 0.23680051381070D-09,
     1         a21= 0.71369114266351D-13,
     2         a22= 0.25045010666356D-09,
     3         a23= 0.72784546444320D-10,
     4         a24= 0.16017135514411D-16,
     5         a25= 0.56562757086698D-10,
     6         a26= 0.28443854062250D-12,
     7         a27= 0.38920900760265D-13,
     8         a28= 0.40317346616717D-21,
     9         a29=-0.92975593753127D-23,
     &         a30=-0.21323943802988D-25,
     1         a31= 0.10007510864939D-25,
     2         a32=-0.15777067572480D-26,
     3         a33= 0.83571296309235D-28)
c
      targ = 1386.d0/t-1.222d0
      targ2 = targ*targ
      targ4 = targ2*targ2
      targ8 = targ4*targ4
      targinv = 1.d0/targ
      targinv2 = targinv*targinv
      targinv3 = targinv2*targinv
      targinv4 = targinv3*targinv
      targinv7 = targinv4*targinv3
      targinv8 = targinv7*targinv
C
      parg = 7.1d0-p*0.60496067755596D-01
      parg2 = parg*parg
      parg3 = parg2*parg
      parg6 = parg3*parg3
c
      pt = parg*targinv
C
      vpt1n=((((((((((((((pt*a33 + a32)*pt + a31)*pt + a30)*
     &      parg6*targinv7 + a29)*parg2*targinv2 + a28)*
     &      parg6*parg6*parg*targinv8*targinv8*targinv2 + a26)*
     &      targinv*targinv4 + a27)*parg3*targinv + targinv3*a25)*
     &      parg + a22)*targinv3 + a23)*targinv2 + targ8*targ2*a24)*
     &      parg + targ4*targ2*a21 + a20 + targinv4*a19)*parg +
     &      (targ8*targ8*a18 + targ2*a17 + a16)*targ + a15 +
     &      targinv3*a14)*parg + (targ2*a13 + a12)*targ + a11 +
     &      (targinv2*a8 + a9)*targinv7 + targinv*a10)*t
C
      END
c
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is pressure in MPa.
C  Second argument is temperature in K.
C  Returns specific enthalpy in kJ/kg
C
      FUNCTION hpt1n(P,T)
C
      IMPLICIT double precision (A-H,O-Z)
C
      parameter(a01=-0.18720692775139D+03,
     2          a02= 0.54083364671138D+03,
     4          a04= 0.21656306556572D+04,
     5          a05=-0.12255145525730D+04,
     6          a06= 0.30266937911228D+03,
     7          a07=-0.42516429081128D+02,
     8          a08= 0.25975485679232D+01,
     9          a09=-0.16303507737913D+01,
     &          a10= 0.27182613947704D+01,
     2          a11= 0.12147472571259D+02,
     3          a13=-0.13971601220485D+02,
     4          a14=-0.10139813560979D+00,
     6          a15= 0.90547896843533D+00,
     7          a17= 0.30487803863262D-01,
     8          a18=-0.84709309503347D-02,
     9          a19=-0.79051996435264D-11,
     1          a20= 0.81058711826910D-01,
     2          a22=-0.32702156038568D-05,
     3          a23= 0.71724465059049D-02,
     4          a24= 0.83376808703816D-03,
     5          a25=-0.91740466143441D-09,
     6          a26= 0.20734169140086D-02,
     7          a27= 0.89603964175208D-05,
     8          a28= 0.66877530790508D-06,
     9          a29= 0.12755771455453D-13,
     &          a30=-0.28710377452428D-15,
     1          a31=-0.64016099916299D-18,
     2          a32= 0.29806124175367D-18,
     3          a33=-0.46640228230284D-19,
     4          a34= 0.24531669269267D-20)
c
      targ = 1386.d0/t-1.222d0
      targ2 = targ*targ
      targ4 = targ2*targ2
      targ8 = targ4*targ4
      targinv = 1.d0/targ
      targinv2 = targinv*targinv
      targinv4 = targinv2*targinv2
      targinv8 = targinv4*targinv4
      targinv9 = targinv8*targinv
C
      parg = 7.1d0-p*0.60496067755596D-01
      parg2 = parg*parg
      parg4 = parg2*parg2
      parg8 = parg4*parg4
C
      pt = parg*targinv
c
      hpt1n = ((((((((((((pt*a34 + a33)*pt + a32)*pt + a31)*
     &   parg4*parg2*targinv9 + targinv2*a30)*parg2 +
     &   a29)*parg8*parg4*parg*targinv8*targinv8*targinv8 +
     &   targinv*a28)*targinv2 + targinv8*a27)*parg2*parg*
     &   targinv4 + targinv9*a26)*parg + targinv2*(targinv*
     &   a24 + targinv4*a23) + targ8*targ*a25)*parg + targ4*targ*
     &   a22 + targinv4*targinv*a20)*parg + targ8*targ8*a19 +
     &   targ2*a18 + targinv4*a15 + a17)*parg + targinv8*
     &   (targinv2*a09 + a10) + targ2*a14 + targinv2*a11 + a13)*parg +
     &   (targ2*a08 + targ*a07 + a06)*targ2 +
     &   targinv2*(targinv*a01 + a02) + targ*a05 + a04
C
      END
c
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is pressure in MPa.
C  Second argument is temperature in K.
C  Returns specific entropy in kJ/(kg*K)
C
      FUNCTION SPT1N(P,T)
C
      IMPLICIT double precision (A-H,O-Z)
C
      PARAMETER (a01=  0.146 329 712 131 67 D0,
     1           a02= -0.845 481 871 691 14 D0,
     2           a03= -0.375 636 036 720 40 D01,
     3           a04=  0.338 551 691 683 85 D01,
     4           a05= -0.957 919 633 878 72 D0,
     5           a06=  0.157 720 385 132 28 D0,
     6           a07= -0.166 164 171 995 01 D-01,
     7           a08=  0.812 146 299 835 68 D-03,
     8           a09=  0.283 190 801 238 04 D-03,
     9           a10= -0.607 063 015 658 74 D-03,
     1           a11= -0.189 900 682 184 19 D-01,
     2           a12= -0.325 297 487 705 05 D-01,
     3           a13= -0.218 417 171 754 14 D-01,
     4           a14= -0.528 383 579 699 30 D-04,
     5           a15= -0.471 843 210 732 67 D-03,
     6           a16= -0.300 017 807 930 26 D-03,
     7           a17=  0.476 613 939 069 87 D-04,
     8           a18= -0.441 418 453 308 46 D-05,
     9           a19= -0.726 949 962 975 94 D-15,
     1           a20= -0.316 796 448 450 54 D-04,
     2           a21= -0.282 707 979 853 12 D-05,
     3           a22= -0.852 051 281 201 03 D-09,
     4           a23= -0.224 252 819 080 00 D-05,
     5           a24= -0.651 712 228 956 01 D-06,
     6           a25= -0.143 417 299 379 24 D-12,
     7           a26= -0.405 169 968 601 17 D-06,
     8           a27= -0.127 343 017 416 41 D-08,
     9           a28= -0.174 248 712 306 34 D-09,
     1           a29= -0.687 621 312 955 31 D-18,
     2           a30=  0.144 783 078 285 21 D-19,
     3           a31=  0.263 357 816 627 95 D-22,
     4           a32= -0.119 476 226 400 71 D-22,
     5           a33=  0.182 280 945 814 04 D-23,
     6           a34= -0.935 370 872 924 58 D-25)
C
      PARAMETER (b01= -0.146 329 712 131 67 D0*2.D0,
     1           b02=  0.845 481 871 691 14 D0,
     3           b04=  0.338 551 691 683 85 D01,
     4           b05= -0.957 919 633 878 72 D0*2.D0,
     5           b06=  0.157 720 385 132 28 D0*3.D0,
     6           b07= -0.166 164 171 995 01 D-01*4.D0,
     7           b08=  0.812 146 299 835 68 D-03*5.D0,
     8           b09= -0.283 190 801 238 04 D-03*9.D0,
     9           b10=  0.607 063 015 658 74 D-03*7.D0,
     1           b11=  0.189 900 682 184 19 D-01,
     3           b13= -0.218 417 171 754 14 D-01,
     4           b14= -0.528 383 579 699 30 D-04*3.D0,
     5           b15=  0.471 843 210 732 67 D-03*3.D0,
     7           b17=  0.476 613 939 069 87 D-04,
     8           b18= -0.441 418 453 308 46 D-05*3.D0,
     9           b19= -0.726 949 962 975 94 D-15*17.D0,
     1           b20=  0.316 796 448 450 54 D-04*4.D0,
     3           b22= -0.852 051 281 201 03 D-09*6.D0,
     4           b23=  0.224 252 819 080 00 D-05*5.D0,
     5           b24=  0.651 712 228 956 01 D-06*2.D0,
     6           b25= -0.143 417 299 379 24 D-12*10.D0,
     7           b26=  0.405 169 968 601 17 D-06*8.D0,
     8           b27=  0.127 343 017 416 41 D-08*11.D0,
     9           b28=  0.174 248 712 306 34 D-09*6.D0,
     1           b29=  0.687 621 312 955 31 D-18*29.D0,
     2           b30= -0.144 783 078 285 21 D-19*31.D0,
     3           b31= -0.263 357 816 627 95 D-22*38.D0,
     4           b32=  0.119 476 226 400 71 D-22*39.D0,
     5           b33= -0.182 280 945 814 04 D-23*40.D0,
     6           b34=  0.935 370 872 924 58 D-25*41.D0)
C
      PARAMETER (tn= 1386.D0,
     1           pn= 16.53D0,
     2           pnq= 1.D0/pn,
     3           r= 0.461526D0)
C
      tau= tn/t
      taun= tau - 1.222d0
      taun2= taun*taun
      taun3= taun2*taun
      taun6= taun3*taun3
      taun12= taun6*taun6
      tauninv= 1.d0/taun
      tauninv2= tauninv*tauninv
      tauninv3= tauninv2*tauninv
      tauninv7= tauninv3*tauninv3*tauninv
C
      pin= 7.1d0 - p*pnq
      pin2= pin*pin
      pin3= pin2*pin
      pin6= pin3*pin3
C
      ga1= tauninv2*(a01 + taun*(a02 + taun*(a03 + taun*(a04 +
     1     taun*(a05 + taun*(a06 + taun*(a07 + taun*a08)))))) +
     2     pin*  (tauninv7*(a09 + taun2*(a10 + taun6*(a11 +
     3           taun*(a12 + taun*(a13 + taun2*a14))))))) +
     4     pin2* tauninv3*(a15 + taun3*(a16 + taun*(a17 +
     5           taun2*(a18 + taun12*taun2*a19))) +
     6     pin*  tauninv*(a20 + taun3*taun*(a21 + taun6*a22) +
     7     pin*  tauninv*(a23 + taun3*(a24 + taun12*a25) +
     8     pin*  tauninv3*(a26 +
     9     pin3* (tauninv3*(a27 + taun3*taun2*a28) +
     1     pin6*pin6*pin*tauninv7*tauninv7*tauninv7*(a29 +
     2     pin2* tauninv2*(a30 +
     3     pin6* tauninv7*(a31 +
     4     pin*  tauninv*(a32 +
     5     pin*  tauninv*(a33 +
     6     pin*  tauninv*a34))) ))) ))) )
C
      gat1= tauninv3*(b01 + taun*(b02 + taun2*(b04 +
     1     taun*(b05 + taun*(b06 + taun*(b07 + taun*b08))))) +
     2     pin*  (tauninv7*(b09 + taun2*(b10 + taun6*(b11 +
     3           taun2*(b13 + taun2*b14)))) +
     4     pin*  tauninv*(b15 + taun3*taun*(b17 +
     5           taun2*(b18 + taun12*taun2*b19)) +
     6     pin*  tauninv*(b20 + taun6*taun3*taun*b22 +
     7     pin*  tauninv*(b23 + taun3*(b24 + taun12*b25) +
     8     pin*  tauninv3*(b26 +
     9     pin3*(tauninv3*(b27 + taun3*taun2*b28) +
     1     pin6*pin6*pin*tauninv7*tauninv7*tauninv7*(b29 +
     2     pin2* tauninv2*(b30 +
     3     pin6* tauninv7*(b31 +
     4     pin*  tauninv*(b32 +
     5     pin*  tauninv*(b33 +
     6     pin*  tauninv*b34))) ))) ))) )))
C
      spt1n= (tau*gat1 - ga1)*r
C
      END
C
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is pressure in MPa.
C  Second argument is temperature in K.
C  Returns speed of sound in m/s
C
      FUNCTION WPT1N(P,T)
C
      IMPLICIT double precision (A-H,O-Z)
C
      PARAMETER (e01=  0.146 329 712 131 67 D0*6.D0,
     1           e02= -0.845 481 871 691 14 D0*2.D0,
     4           e05= -0.957 919 633 878 72 D0*2.D0,
     5           e06=  0.157 720 385 132 28 D0*6.D0,
     6           e07= -0.166 164 171 995 01 D-01*12.D0,
     7           e08=  0.812 146 299 835 68 D-03*20.D0,
     8           e09=  0.283 190 801 238 04 D-03*90.D0,
     9           e10= -0.607 063 015 658 74 D-03*56.D0,
     1           e11= -0.189 900 682 184 19 D-01*2.D0,
     4           e14= -0.528 383 579 699 30 D-04*6.D0,
     5           e15= -0.471 843 210 732 67 D-03*12.D0,
     8           e18= -0.441 418 453 308 46 D-05*6.D0,
     9           e19= -0.726 949 962 975 94 D-15*272.D0,
     1           e20= -0.316 796 448 450 54 D-04*20.D0,
     3           e22= -0.852 051 281 201 03 D-09*30.D0,
     4           e23= -0.224 252 819 080 00 D-05*30.D0,
     5           e24= -0.651 712 228 956 01 D-06*6.D0,
     6           e25= -0.143 417 299 379 24 D-12*90.D0,
     7           e26= -0.405 169 968 601 17 D-06*72.D0,
     8           e27= -0.127 343 017 416 41 D-08*132.D0,
     9           e28= -0.174 248 712 306 34 D-09*42.D0,
     1           e29= -0.687 621 312 955 31 D-18*870.D0,
     2           e30=  0.144 783 078 285 21 D-19*992.D0,
     3           e31=  0.263 357 816 627 95 D-22*1482.D0,
     4           e32= -0.119 476 226 400 71 D-22*1560.D0,
     5           e33=  0.182 280 945 814 04 D-23*1640.D0,
     6           e34= -0.935 370 872 924 58 D-25*1722.D0)
C
      PARAMETER (c15= -0.471 843 210 732 67 D-03*2.D0,
     6           c16= -0.300 017 807 930 26 D-03*2.D0,
     7           c17=  0.476 613 939 069 87 D-04*2.D0,
     8           c18= -0.441 418 453 308 46 D-05*2.D0,
     9           c19= -0.726 949 962 975 94 D-15*2.D0,
     1           c20= -0.316 796 448 450 54 D-04*6.D0,
     2           c21= -0.282 707 979 853 12 D-05*6.D0,
     3           c22= -0.852 051 281 201 03 D-09*6.D0,
     4           c23= -0.224 252 819 080 00 D-05*12.D0,
     5           c24= -0.651 712 228 956 01 D-06*12.D0,
     6           c25= -0.143 417 299 379 24 D-12*12.D0,
     7           c26= -0.405 169 968 601 17 D-06*20.D0,
     8           c27= -0.127 343 017 416 41 D-08*56.D0,
     9           c28= -0.174 248 712 306 34 D-09*56.D0,
     1           c29= -0.687 621 312 955 31 D-18*420.D0,
     2           c30=  0.144 783 078 285 21 D-19*506.D0,
     3           c31=  0.263 357 816 627 95 D-22*812.D0,
     4           c32= -0.119 476 226 400 71 D-22*870.D0,
     5           c33=  0.182 280 945 814 04 D-23*930.D0,
     6           c34= -0.935 370 872 924 58 D-25*992.D0)
C
      PARAMETER (b09=  0.283 190 801 238 04 D-03,
     9           b10= -0.607 063 015 658 74 D-03,
     1           b11= -0.189 900 682 184 19 D-01,
     2           b12= -0.325 297 487 705 05 D-01,
     3           b13= -0.218 417 171 754 14 D-01,
     4           b14= -0.528 383 579 699 30 D-04,
     5           b15= -0.471 843 210 732 67 D-03*2.D0,
     6           b16= -0.300 017 807 930 26 D-03*2.D0,
     7           b17=  0.476 613 939 069 87 D-04*2.D0,
     8           b18= -0.441 418 453 308 46 D-05*2.D0,
     9           b19= -0.726 949 962 975 94 D-15*2.D0,
     1           b20= -0.316 796 448 450 54 D-04*3.D0,
     2           b21= -0.282 707 979 853 12 D-05*3.D0,
     3           b22= -0.852 051 281 201 03 D-09*3.D0,
     4           b23= -0.224 252 819 080 00 D-05*4.D0,
     5           b24= -0.651 712 228 956 01 D-06*4.D0,
     6           b25= -0.143 417 299 379 24 D-12*4.D0,
     7           b26= -0.405 169 968 601 17 D-06*5.D0,
     8           b27= -0.127 343 017 416 41 D-08*8.D0,
     9           b28= -0.174 248 712 306 34 D-09*8.D0,
     1           b29= -0.687 621 312 955 31 D-18*21.D0,
     2           b30=  0.144 783 078 285 21 D-19*23.D0,
     3           b31=  0.263 357 816 627 95 D-22*29.D0,
     4           b32= -0.119 476 226 400 71 D-22*30.D0,
     5           b33=  0.182 280 945 814 04 D-23*31.D0,
     6           b34= -0.935 370 872 924 58 D-25*32.D0)
C
      PARAMETER (f09= -0.283 190 801 238 04 D-03*9.D0,
     9           f10=  0.607 063 015 658 74 D-03*7.D0,
     1           f11=  0.189 900 682 184 19 D-01*1.D0,
     3           f13= -0.218 417 171 754 14 D-01,
     4           f14= -0.528 383 579 699 30 D-04*3.D0,
     5           f15=  0.471 843 210 732 67 D-03*6.D0,
     7           f17=  0.476 613 939 069 87 D-04*2.D0,
     8           f18= -0.441 418 453 308 46 D-05*6.D0,
     9           f19= -0.726 949 962 975 94 D-15*34.D0,
     1           f20=  0.316 796 448 450 54 D-04*12.D0,
     3           f22= -0.852 051 281 201 03 D-09*18.D0,
     4           f23=  0.224 252 819 080 00 D-05*20.D0,
     5           f24=  0.651 712 228 956 01 D-06*8.D0,
     6           f25= -0.143 417 299 379 24 D-12*40.D0,
     7           f26=  0.405 169 968 601 17 D-06*40.D0,
     8           f27=  0.127 343 017 416 41 D-08*88.D0,
     9           f28=  0.174 248 712 306 34 D-09*48.D0,
     1           f29=  0.687 621 312 955 31 D-18*609.D0,
     2           f30= -0.144 783 078 285 21 D-19*713.D0,
     3           f31= -0.263 357 816 627 95 D-22*1102.D0,
     4           f32=  0.119 476 226 400 71 D-22*1170.D0,
     5           f33= -0.182 280 945 814 04 D-23*1240.D0,
     6           f34=  0.935 370 872 924 58 D-25*1312.D0)
C
      PARAMETER (tn= 1386.D0,
     1           pn= 16.53D0,
     2           pnq= 1.D0/pn,
     3           r= 0.461526D0,
     4           rn= r*1000.0D0)
C
      tau= tn/t
      taun= tau - 1.222d0
      taun2= taun*taun
      taun3= taun2*taun
      taun6= taun3*taun3
      taun12= taun6*taun6
      tauninv= 1.d0/taun
      tauninv2= tauninv*tauninv
      tauninv3= tauninv2*tauninv
      tauninv7= tauninv3*tauninv3*tauninv
C
      pin= 7.1d0 - p*pnq
      pin2= pin*pin
      pin3= pin2*pin
      pin6= pin3*pin3
C
      gatt1= tauninv2*(tauninv2*(e01 + taun*(e02 + taun3*(e05 +
     1       taun*(e06 + taun*(e07 + taun*e08)))) +
     2       pin*  (tauninv7*(e09 + taun2*(e10 + taun6*(e11 +
     3             taun2*taun2*e14))))) +
     4       pin2* tauninv3*(e15 + taun6*(e18 +
     5             taun12*taun2*e19) +
     6       pin*  tauninv*(e20 + taun6*taun3*taun*e22 +
     7       pin*  tauninv*(e23 + taun3*(e24 + taun12*e25) +
     8       pin*  tauninv3*(e26 +
     9       pin3*(tauninv3*(e27 + taun3*taun2*e28) +
     1       pin6*pin6*pin*tauninv7*tauninv7*tauninv7*(e29 +
     2       pin2* tauninv2*(e30 +
     3       pin6* tauninv7*(e31 +
     4       pin*  tauninv*(e32 +
     5       pin*  tauninv*(e33 +
     6       pin*  tauninv*e34))) ))) ))) ))
C
      gapp1= tauninv3*(c15 + taun3*(c16 + taun*(c17 +
     5       taun2*(c18 + taun12*taun2*c19))) +
     6       pin*  tauninv*(c20 + taun3*taun*(c21 + taun6*c22) +
     7       pin*  tauninv*(c23 + taun3*(c24 + taun12*c25) +
     8       pin*  tauninv3*(c26 +
     9       pin3* (tauninv3*(c27 + taun3*taun2*c28) +
     1       pin6*pin6*pin*tauninv7*tauninv7*tauninv7*(c29 +
     2       pin2* tauninv2*(c30 +
     3       pin6* tauninv7*(c31 +
     4       pin*  tauninv*(c32 +
     5       pin*  tauninv*(c33 +
     6       pin*  tauninv*c34))) ))) ))) )
C
      gap1= -(tauninv7*(tauninv2*b09 + b10 + taun6*(b11 +
     1      taun*(b12 + taun*(b13 + taun2*b14)))) +
     2      pin*  tauninv3*(b15 + taun3*(b16 + taun*(b17 +
     3            taun2*(b18 + taun12*taun2*b19))) +
     4      pin*  tauninv*(b20 + taun3*taun*(b21 + taun6*b22) +
     5      pin*  tauninv*(b23 + taun3*(b24 + taun12*b25) +
     6      pin*  tauninv3*(b26 +
     7      pin3* (tauninv3*(b27 + taun3*taun2*b28) +
     8      pin6*pin6*pin*tauninv7*tauninv7*tauninv7*(b29 +
     9      pin2* tauninv2*(b30 +
     1      pin6* tauninv7*(b31 +
     2      pin*  tauninv*(b32 +
     3      pin*  tauninv*(b33 +
     4      pin*  tauninv*b34))) ))) ))) ))
C
      gapt1= -(tauninv*(tauninv7*(tauninv2*f09 + f10 + taun6*(f11 +
     3       taun2*(f13 + taun2*f14))) +
     4       pin*  tauninv3*(f15 + taun3*taun*(f17 +
     5             taun2*(f18 + taun12*taun2*f19)) +
     6       pin*  tauninv*(f20 + taun6*taun3*taun*f22 +
     7       pin*  tauninv*(f23 + taun3*(f24 + taun12*f25) +
     8       pin*  tauninv3*(f26 +
     9       pin3*(tauninv3*(f27 + taun3*taun2*f28) +
     1       pin6*pin6*pin*tauninv7*tauninv7*tauninv7*(f29 +
     2       pin2* tauninv2*(f30 +
     3       pin6* tauninv7*(f31 +
     4       pin*  tauninv*(f32 +
     5       pin*  tauninv*(f33 +
     6       pin*  tauninv*f34))) ))) ))) )))
C
      wpt1n= DSQRT(rn*t*(gap1*gap1)/(((gap1 - tau*gapt1)*
     1       (gap1 - tau*gapt1))/(tau*tau*gatt1) - gapp1))
C
      END
C
c
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is pressure in MPa.
C  Second argument is temperature in K.
C  Returns specific isobaric heat capacity in kJ/(kg*K)
C
      FUNCTION CPPT1N(P,T)
C
      IMPLICIT double precision (A-H,O-Z)
C
      PARAMETER (e01=  0.146 329 712 131 67 D0*6.D0,
     1           e02= -0.845 481 871 691 14 D0*2.D0,
     4           e05= -0.957 919 633 878 72 D0*2.D0,
     5           e06=  0.157 720 385 132 28 D0*6.D0,
     6           e07= -0.166 164 171 995 01 D-01*12.D0,
     7           e08=  0.812 146 299 835 68 D-03*20.D0,
     8           e09=  0.283 190 801 238 04 D-03*90.D0,
     9           e10= -0.607 063 015 658 74 D-03*56.D0,
     1           e11= -0.189 900 682 184 19 D-01*2.D0,
     4           e14= -0.528 383 579 699 30 D-04*6.D0,
     5           e15= -0.471 843 210 732 67 D-03*12.D0,
     8           e18= -0.441 418 453 308 46 D-05*6.D0,
     9           e19= -0.726 949 962 975 94 D-15*272.D0,
     1           e20= -0.316 796 448 450 54 D-04*20.D0,
     3           e22= -0.852 051 281 201 03 D-09*30.D0,
     4           e23= -0.224 252 819 080 00 D-05*30.D0,
     5           e24= -0.651 712 228 956 01 D-06*6.D0,
     6           e25= -0.143 417 299 379 24 D-12*90.D0,
     7           e26= -0.405 169 968 601 17 D-06*72.D0,
     8           e27= -0.127 343 017 416 41 D-08*132.D0,
     9           e28= -0.174 248 712 306 34 D-09*42.D0,
     1           e29= -0.687 621 312 955 31 D-18*870.D0,
     2           e30=  0.144 783 078 285 21 D-19*992.D0,
     3           e31=  0.263 357 816 627 95 D-22*1482.D0,
     4           e32= -0.119 476 226 400 71 D-22*1560.D0,
     5           e33=  0.182 280 945 814 04 D-23*1640.D0,
     6           e34= -0.935 370 872 924 58 D-25*1722.D0)
C
      PARAMETER (tn= 1386.0D0,
     1           pn= 16.53D0,
     2           pnq= 1.D0/pn,
     3           r= 0.461526D0)
C
      tau= tn/t
      taun= tau - 1.222D0
      taun2= taun*taun
      taun3= taun2*taun
      taun6= taun3*taun3
      taun12= taun6*taun6
      tauninv= 1.d0/taun
      tauninv2= tauninv*tauninv
      tauninv3= tauninv2*tauninv
      tauninv7= tauninv3*tauninv3*tauninv
C
      pin= 7.1d0 - p*pnq
      pin2= pin*pin
      pin3= pin2*pin
      pin6= pin3*pin3
C
      cppt1n= (tauninv2*(tauninv2*(e01 + taun*(e02 + taun3*(e05 +
     1       taun*(e06 + taun*(e07 + taun*e08)))) +
     2       pin*  (tauninv7*(e09 + taun2*(e10 + taun6*(e11 +
     3             taun2*taun2*e14))))) +
     4       pin2* tauninv3*(e15 + taun6*(e18 +
     5             taun12*taun2*e19) +
     6       pin*  tauninv*(e20 + taun6*taun3*taun*e22 +
     7       pin*  tauninv*(e23 + taun3*(e24 + taun12*e25) +
     8       pin*  tauninv3*(e26 +
     9       pin3*(tauninv3*(e27 + taun3*taun2*e28) +
     1       pin6*pin6*pin*tauninv7*tauninv7*tauninv7*(e29 +
     2       pin2* tauninv2*(e30 +
     3       pin6* tauninv7*(e31 +
     4       pin*  tauninv*(e32 +
     5       pin*  tauninv*(e33 +
     6       pin*  tauninv*e34))))))))))))*(-tau*tau*r)
C
      END
C
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is pressure in MPa.
C  Second argument is temperature in K.
C  Returns isochoric heat capacity in kJ/(kg * K)
C
      FUNCTION CVPT1N(P,T)
C
      IMPLICIT double precision (A-H,O-Z)
C
      PARAMETER (e01=  0.146 329 712 131 67 D0*6.D0,
     1           e02= -0.845 481 871 691 14 D0*2.D0,
     4           e05= -0.957 919 633 878 72 D0*2.D0,
     5           e06=  0.157 720 385 132 28 D0*6.D0,
     6           e07= -0.166 164 171 995 01 D-01*12.D0,
     7           e08=  0.812 146 299 835 68 D-03*20.D0,
     8           e09=  0.283 190 801 238 04 D-03*90.D0,
     9           e10= -0.607 063 015 658 74 D-03*56.D0,
     1           e11= -0.189 900 682 184 19 D-01*2.D0,
     4           e14= -0.528 383 579 699 30 D-04*6.D0,
     5           e15= -0.471 843 210 732 67 D-03*12.D0,
     8           e18= -0.441 418 453 308 46 D-05*6.D0,
     9           e19= -0.726 949 962 975 94 D-15*272.D0,
     1           e20= -0.316 796 448 450 54 D-04*20.D0,
     3           e22= -0.852 051 281 201 03 D-09*30.D0,
     4           e23= -0.224 252 819 080 00 D-05*30.D0,
     5           e24= -0.651 712 228 956 01 D-06*6.D0,
     6           e25= -0.143 417 299 379 24 D-12*90.D0,
     7           e26= -0.405 169 968 601 17 D-06*72.D0,
     8           e27= -0.127 343 017 416 41 D-08*132.D0,
     9           e28= -0.174 248 712 306 34 D-09*42.D0,
     1           e29= -0.687 621 312 955 31 D-18*870.D0,
     2           e30=  0.144 783 078 285 21 D-19*992.D0,
     3           e31=  0.263 357 816 627 95 D-22*1482.D0,
     4           e32= -0.119 476 226 400 71 D-22*1560.D0,
     5           e33=  0.182 280 945 814 04 D-23*1640.D0,
     6           e34= -0.935 370 872 924 58 D-25*1722.D0)
C
      PARAMETER (c15= -0.471 843 210 732 67 D-03*2.D0,
     6           c16= -0.300 017 807 930 26 D-03*2.D0,
     7           c17=  0.476 613 939 069 87 D-04*2.D0,
     8           c18= -0.441 418 453 308 46 D-05*2.D0,
     9           c19= -0.726 949 962 975 94 D-15*2.D0,
     1           c20= -0.316 796 448 450 54 D-04*6.D0,
     2           c21= -0.282 707 979 853 12 D-05*6.D0,
     3           c22= -0.852 051 281 201 03 D-09*6.D0,
     4           c23= -0.224 252 819 080 00 D-05*12.D0,
     5           c24= -0.651 712 228 956 01 D-06*12.D0,
     6           c25= -0.143 417 299 379 24 D-12*12.D0,
     7           c26= -0.405 169 968 601 17 D-06*20.D0,
     8           c27= -0.127 343 017 416 41 D-08*56.D0,
     9           c28= -0.174 248 712 306 34 D-09*56.D0,
     1           c29= -0.687 621 312 955 31 D-18*420.D0,
     2           c30=  0.144 783 078 285 21 D-19*506.D0,
     3           c31=  0.263 357 816 627 95 D-22*812.D0,
     4           c32= -0.119 476 226 400 71 D-22*870.D0,
     5           c33=  0.182 280 945 814 04 D-23*930.D0,
     6           c34= -0.935 370 872 924 58 D-25*992.D0)
C
      PARAMETER (b09=  0.283 190 801 238 04 D-03,
     9           b10= -0.607 063 015 658 74 D-03,
     1           b11= -0.189 900 682 184 19 D-01,
     2           b12= -0.325 297 487 705 05 D-01,
     3           b13= -0.218 417 171 754 14 D-01,
     4           b14= -0.528 383 579 699 30 D-04,
     5           b15= -0.471 843 210 732 67 D-03*2.D0,
     6           b16= -0.300 017 807 930 26 D-03*2.D0,
     7           b17=  0.476 613 939 069 87 D-04*2.D0,
     8           b18= -0.441 418 453 308 46 D-05*2.D0,
     9           b19= -0.726 949 962 975 94 D-15*2.D0,
     1           b20= -0.316 796 448 450 54 D-04*3.D0,
     2           b21= -0.282 707 979 853 12 D-05*3.D0,
     3           b22= -0.852 051 281 201 03 D-09*3.D0,
     4           b23= -0.224 252 819 080 00 D-05*4.D0,
     5           b24= -0.651 712 228 956 01 D-06*4.D0,
     6           b25= -0.143 417 299 379 24 D-12*4.D0,
     7           b26= -0.405 169 968 601 17 D-06*5.D0,
     8           b27= -0.127 343 017 416 41 D-08*8.D0,
     9           b28= -0.174 248 712 306 34 D-09*8.D0,
     1           b29= -0.687 621 312 955 31 D-18*21.D0,
     2           b30=  0.144 783 078 285 21 D-19*23.D0,
     3           b31=  0.263 357 816 627 95 D-22*29.D0,
     4           b32= -0.119 476 226 400 71 D-22*30.D0,
     5           b33=  0.182 280 945 814 04 D-23*31.D0,
     6           b34= -0.935 370 872 924 58 D-25*32.D0)
C
      PARAMETER (f09= -0.283 190 801 238 04 D-03*9.D0,
     9           f10=  0.607 063 015 658 74 D-03*7.D0,
     1           f11=  0.189 900 682 184 19 D-01*1.D0,
     3           f13= -0.218 417 171 754 14 D-01,
     4           f14= -0.528 383 579 699 30 D-04*3.D0,
     5           f15=  0.471 843 210 732 67 D-03*6.D0,
     7           f17=  0.476 613 939 069 87 D-04*2.D0,
     8           f18= -0.441 418 453 308 46 D-05*6.D0,
     9           f19= -0.726 949 962 975 94 D-15*34.D0,
     1           f20=  0.316 796 448 450 54 D-04*12.D0,
     3           f22= -0.852 051 281 201 03 D-09*18.D0,
     4           f23=  0.224 252 819 080 00 D-05*20.D0,
     5           f24=  0.651 712 228 956 01 D-06*8.D0,
     6           f25= -0.143 417 299 379 24 D-12*40.D0,
     7           f26=  0.405 169 968 601 17 D-06*40.D0,
     8           f27=  0.127 343 017 416 41 D-08*88.D0,
     9           f28=  0.174 248 712 306 34 D-09*48.D0,
     1           f29=  0.687 621 312 955 31 D-18*609.D0,
     2           f30= -0.144 783 078 285 21 D-19*713.D0,
     3           f31= -0.263 357 816 627 95 D-22*1102.D0,
     4           f32=  0.119 476 226 400 71 D-22*1170.D0,
     5           f33= -0.182 280 945 814 04 D-23*1240.D0,
     6           f34=  0.935 370 872 924 58 D-25*1312.D0)
C
      PARAMETER (tn= 1386.D0,
     1           pn= 16.53D0,
     2           pnq= 1.D0/pn,
     3           r= 0.461526D0)
C
      tau= tn/t
      taun= tau - 1.222d0
      taun2= taun*taun
      taun3= taun2*taun
      taun6= taun3*taun3
      taun12= taun6*taun6
      tauninv= 1.d0/taun
      tauninv2= tauninv*tauninv
      tauninv3= tauninv2*tauninv
      tauninv7= tauninv3*tauninv3*tauninv
C
      pin= 7.1d0 - p*pnq
      pin2= pin*pin
      pin3= pin2*pin
      pin6= pin3*pin3
C
      gatt1= tauninv2*(tauninv2*(e01 + taun*(e02 + taun3*(e05 +
     1       taun*(e06 + taun*(e07 + taun*e08)))) +
     2       pin*  (tauninv7*(e09 + taun2*(e10 + taun6*(e11 +
     3             taun2*taun2*e14))))) +
     4       pin2* tauninv3*(e15 + taun6*(e18 +
     5             taun12*taun2*e19) +
     6       pin*  tauninv*(e20 + taun6*taun3*taun*e22 +
     7       pin*  tauninv*(e23 + taun3*(e24 + taun12*e25) +
     8       pin*  tauninv3*(e26 +
     9       pin3*(tauninv3*(e27 + taun3*taun2*e28) +
     1       pin6*pin6*pin*tauninv7*tauninv7*tauninv7*(e29 +
     2       pin2* tauninv2*(e30 +
     3       pin6* tauninv7*(e31 +
     4       pin*  tauninv*(e32 +
     5       pin*  tauninv*(e33 +
     6       pin*  tauninv*e34))) ))) ))) ))
C
      gapp1= tauninv3*(c15 + taun3*(c16 + taun*(c17 +
     5       taun2*(c18 + taun12*taun2*c19))) +
     6       pin*  tauninv*(c20 + taun3*taun*(c21 + taun6*c22) +
     7       pin*  tauninv*(c23 + taun3*(c24 + taun12*c25) +
     8       pin*  tauninv3*(c26 +
     9       pin3* (tauninv3*(c27 + taun3*taun2*c28) +
     1       pin6*pin6*pin*tauninv7*tauninv7*tauninv7*(c29 +
     2       pin2* tauninv2*(c30 +
     3       pin6* tauninv7*(c31 +
     4       pin*  tauninv*(c32 +
     5       pin*  tauninv*(c33 +
     6       pin*  tauninv*c34))) ))) ))) )
C
      gap1= -(tauninv7*(tauninv2*b09 + b10 + taun6*(b11 +
     1      taun*(b12 + taun*(b13 + taun2*b14)))) +
     2      pin*  tauninv3*(b15 + taun3*(b16 + taun*(b17 +
     3            taun2*(b18 + taun12*taun2*b19))) +
     4      pin*  tauninv*(b20 + taun3*taun*(b21 + taun6*b22) +
     5      pin*  tauninv*(b23 + taun3*(b24 + taun12*b25) +
     6      pin*  tauninv3*(b26 +
     7      pin3* (tauninv3*(b27 + taun3*taun2*b28) +
     8      pin6*pin6*pin*tauninv7*tauninv7*tauninv7*(b29 +
     9      pin2* tauninv2*(b30 +
     1      pin6* tauninv7*(b31 +
     2      pin*  tauninv*(b32 +
     3      pin*  tauninv*(b33 +
     4      pin*  tauninv*b34))) ))) ))) ))
C
      gapt1= -(tauninv*(tauninv7*(tauninv2*f09 + f10 + taun6*(f11 +
     3       taun2*(f13 + taun2*f14))) +
     4       pin*  tauninv3*(f15 + taun3*taun*(f17 +
     5             taun2*(f18 + taun12*taun2*f19)) +
     6       pin*  tauninv*(f20 + taun6*taun3*taun*f22 +
     7       pin*  tauninv*(f23 + taun3*(f24 + taun12*f25) +
     8       pin*  tauninv3*(f26 +
     9       pin3*(tauninv3*(f27 + taun3*taun2*f28) +
     1       pin6*pin6*pin*tauninv7*tauninv7*tauninv7*(f29 +
     2       pin2* tauninv2*(f30 +
     3       pin6* tauninv7*(f31 +
     4       pin*  tauninv*(f32 +
     5       pin*  tauninv*(f33 +
     6       pin*  tauninv*f34))) ))) ))) )))
C
      cvpt1n= r*(-tau*tau*gatt1+(gap1-tau*gapt1)*(gap1-tau*gapt1)/
     1            gapp1)
C
      END
C
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is pressure in MPa.
C  Second argument is temperature in K.
C  Returns specific Gibbs free energy in kJ/kg
C
      FUNCTION GPT1N(P,T)
C
      IMPLICIT double precision (A-H,O-Z)
C
      PARAMETER (a01=  0.146 329 712 131 67 D0,
     1           a02= -0.845 481 871 691 14 D0,
     2           a03= -0.375 636 036 720 40 D01,
     3           a04=  0.338 551 691 683 85 D01,
     4           a05= -0.957 919 633 878 72 D0,
     5           a06=  0.157 720 385 132 28 D0,
     6           a07= -0.166 164 171 995 01 D-01,
     7           a08=  0.812 146 299 835 68 D-03,
     8           a09=  0.283 190 801 238 04 D-03,
     9           a10= -0.607 063 015 658 74 D-03,
     1           a11= -0.189 900 682 184 19 D-01,
     2           a12= -0.325 297 487 705 05 D-01,
     3           a13= -0.218 417 171 754 14 D-01,
     4           a14= -0.528 383 579 699 30 D-04,
     5           a15= -0.471 843 210 732 67 D-03,
     6           a16= -0.300 017 807 930 26 D-03,
     7           a17=  0.476 613 939 069 87 D-04,
     8           a18= -0.441 418 453 308 46 D-05,
     9           a19= -0.726 949 962 975 94 D-15,
     1           a20= -0.316 796 448 450 54 D-04,
     2           a21= -0.282 707 979 853 12 D-05,
     3           a22= -0.852 051 281 201 03 D-09,
     4           a23= -0.224 252 819 080 00 D-05,
     5           a24= -0.651 712 228 956 01 D-06,
     6           a25= -0.143 417 299 379 24 D-12,
     7           a26= -0.405 169 968 601 17 D-06,
     8           a27= -0.127 343 017 416 41 D-08,
     9           a28= -0.174 248 712 306 34 D-09,
     1           a29= -0.687 621 312 955 31 D-18,
     2           a30=  0.144 783 078 285 21 D-19,
     3           a31=  0.263 357 816 627 95 D-22,
     4           a32= -0.119 476 226 400 71 D-22,
     5           a33=  0.182 280 945 814 04 D-23,
     6           a34= -0.935 370 872 924 58 D-25)
C
      PARAMETER (tn= 1386.D0,
     1           pn= 16.53D0,
     2           pnq= 1.0D0/pn,
     3           r= 0.461526D0)
C
      taun= tn/t - 1.222d0
      taun2= taun*taun
      taun3= taun2*taun
      taun6= taun3*taun3
      taun12= taun6*taun6
      tauninv= 1.d0/taun
      tauninv2= tauninv*tauninv
      tauninv3= tauninv2*tauninv
      tauninv7= tauninv3*tauninv3*tauninv
C
      pin= 7.1d0 - p*pnq
      pin2= pin*pin
      pin3= pin2*pin
      pin6= pin3*pin3
C
      gpt1n= (tauninv2*(a01 + taun*(a02 + taun*(a03 + taun*(a04 +
     1     taun*(a05 + taun*(a06 + taun*(a07 + taun*a08)))))) +
     2     pin*  (tauninv7*(a09 + taun2*(a10 + taun6*(a11 +
     3           taun*(a12 + taun*(a13 + taun2*a14))))))) +
     4     pin2* tauninv3*(a15 + taun3*(a16 + taun*(a17 +
     5           taun2*(a18 + taun12*taun2*a19))) +
     6     pin*  tauninv*(a20 + taun3*taun*(a21 + taun6*a22) +
     7     pin*  tauninv*(a23 + taun3*(a24 + taun12*a25) +
     8     pin*  tauninv3*(a26 +
     9     pin3* (tauninv3*(a27 + taun3*taun2*a28) +
     1     pin6*pin6*pin*tauninv7*tauninv7*tauninv7*(a29 +
     2     pin2* tauninv2*(a30 +
     3     pin6* tauninv7*(a31 +
     4     pin*  tauninv*(a32 +
     5     pin*  tauninv*(a33 +
     6     pin*  tauninv*a34))) ))) ))) ))*r*t
C
      END
C
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is pressure in MPa.
C  Second argument is temperature in K.
C  Returns specific internal energy in kJ/kg
C
      FUNCTION UPT1N(P,T)
C
      IMPLICIT double precision (A-H,O-Z)
C
      PARAMETER (b09=  0.283 190 801 238 04 D-03,
     9           b10= -0.607 063 015 658 74 D-03,
     1           b11= -0.189 900 682 184 19 D-01,
     2           b12= -0.325 297 487 705 05 D-01,
     3           b13= -0.218 417 171 754 14 D-01,
     4           b14= -0.528 383 579 699 30 D-04,
     5           b15= -0.471 843 210 732 67 D-03*2.D0,
     6           b16= -0.300 017 807 930 26 D-03*2.D0,
     7           b17=  0.476 613 939 069 87 D-04*2.D0,
     8           b18= -0.441 418 453 308 46 D-05*2.D0,
     9           b19= -0.726 949 962 975 94 D-15*2.D0,
     1           b20= -0.316 796 448 450 54 D-04*3.D0,
     2           b21= -0.282 707 979 853 12 D-05*3.D0,
     3           b22= -0.852 051 281 201 03 D-09*3.D0,
     4           b23= -0.224 252 819 080 00 D-05*4.D0,
     5           b24= -0.651 712 228 956 01 D-06*4.D0,
     6           b25= -0.143 417 299 379 24 D-12*4.D0,
     7           b26= -0.405 169 968 601 17 D-06*5.D0,
     8           b27= -0.127 343 017 416 41 D-08*8.D0,
     9           b28= -0.174 248 712 306 34 D-09*8.D0,
     1           b29= -0.687 621 312 955 31 D-18*21.D0,
     2           b30=  0.144 783 078 285 21 D-19*23.D0,
     3           b31=  0.263 357 816 627 95 D-22*29.D0,
     4           b32= -0.119 476 226 400 71 D-22*30.D0,
     5           b33=  0.182 280 945 814 04 D-23*31.D0,
     6           b34= -0.935 370 872 924 58 D-25*32.D0)
C
      PARAMETER (d01= -0.146 329 712 131 67 D0*2.D0,
     1           d02=  0.845 481 871 691 14 D0,
     3           d04=  0.338 551 691 683 85 D01,
     4           d05= -0.957 919 633 878 72 D0*2.D0,
     5           d06=  0.157 720 385 132 28 D0*3.D0,
     6           d07= -0.166 164 171 995 01 D-01*4.D0,
     7           d08=  0.812 146 299 835 68 D-03*5.D0,
     8           d09= -0.283 190 801 238 04 D-03*9.D0,
     9           d10=  0.607 063 015 658 74 D-03*7.D0,
     1           d11=  0.189 900 682 184 19 D-01,
     3           d13= -0.218 417 171 754 14 D-01,
     4           d14= -0.528 383 579 699 30 D-04*3.D0,
     5           d15=  0.471 843 210 732 67 D-03*3.D0,
     7           d17=  0.476 613 939 069 87 D-04,
     8           d18= -0.441 418 453 308 46 D-05*3.D0,
     9           d19= -0.726 949 962 975 94 D-15*17.D0,
     1           d20=  0.316 796 448 450 54 D-04*4.D0,
     3           d22= -0.852 051 281 201 03 D-09*6.D0,
     4           d23=  0.224 252 819 080 00 D-05*5.D0,
     5           d24=  0.651 712 228 956 01 D-06*2.D0,
     6           d25= -0.143 417 299 379 24 D-12*10.D0,
     7           d26=  0.405 169 968 601 17 D-06*8.D0,
     8           d27=  0.127 343 017 416 41 D-08*11.D0,
     9           d28=  0.174 248 712 306 34 D-09*6.D0,
     1           d29=  0.687 621 312 955 31 D-18*29.D0,
     2           d30= -0.144 783 078 285 21 D-19*31.D0,
     3           d31= -0.263 357 816 627 95 D-22*38.D0,
     4           d32=  0.119 476 226 400 71 D-22*39.D0,
     5           d33= -0.182 280 945 814 04 D-23*40.D0,
     6           d34=  0.935 370 872 924 58 D-25*41.D0)
C
       PARAMETER (tn= 1386.D0,
     1            pn= 16.53D0,
     2            pnq= 1.D0/pn,
     3            r= 0.461526D0)
C
      tau= tn/t
      taun= tau - 1.222d0
      taun2= taun*taun
      taun3= taun2*taun
      taun6= taun3*taun3
      taun12= taun6*taun6
      tauninv= 1.d0/taun
      tauninv2= tauninv*tauninv
      tauninv3= tauninv2*tauninv
      tauninv7= tauninv3*tauninv3*tauninv
C
      pi= p*pnq
      pin= 7.1d0 - pi
      pin2= pin*pin
      pin3= pin2*pin
      pin6= pin3*pin3
C
      gap1= -(tauninv7*(tauninv2*b09 + b10 + taun6*(b11 +
     1      taun*(b12 + taun*(b13 + taun2*b14)))) +
     2      pin*  tauninv3*(b15 + taun3*(b16 + taun*(b17 +
     3            taun2*(b18 + taun12*taun2*b19))) +
     4      pin*  tauninv*(b20 + taun3*taun*(b21 + taun6*b22) +
     5      pin*  tauninv*(b23 + taun3*(b24 + taun12*b25) +
     6      pin*  tauninv3*(b26 +
     7      pin3* (tauninv3*(b27 + taun3*taun2*b28) +
     8      pin6*pin6*pin*tauninv7*tauninv7*tauninv7*(b29 +
     9      pin2* tauninv2*(b30 +
     1      pin6* tauninv7*(b31 +
     2      pin*  tauninv*(b32 +
     3      pin*  tauninv*(b33 +
     4      pin*  tauninv*b34))) ))) ))) ))
C
      gat1= tauninv3*(d01 + taun*(d02 + taun2*(d04 +
     1     taun*(d05 + taun*(d06 + taun*(d07 + taun*d08))))) +
     2     pin*  (tauninv7*(d09 + taun2*(d10 + taun6*(d11 +
     3           taun2*(d13 + taun2*d14)))) +
     4     pin*  tauninv*(d15 + taun3*taun*(d17 +
     5           taun2*(d18 + taun12*taun2*d19)) +
     6     pin*  tauninv*(d20 + taun6*taun3*taun*d22 +
     7     pin*  tauninv*(d23 + taun3*(d24 + taun12*d25) +
     8     pin*  tauninv3*(d26 +
     9     pin3*(tauninv3*(d27 + taun3*taun2*d28) +
     1     pin6*pin6*pin*tauninv7*tauninv7*tauninv7*(d29 +
     2     pin2* tauninv2*(d30 +
     3     pin6* tauninv7*(d31 +
     4     pin*  tauninv*(d32 +
     5     pin*  tauninv*(d33 +
     6     pin*  tauninv*d34))) ))) ))) )))
C
      upt1n= r*t*(tau*gat1 - pi*gap1)
C
      END
C
c
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is pressure in MPa.
C  Second argument is temperature in K.
C  Returns speed of sound in m/s
C
      FUNCTION WPT2N(P,T)
C
      IMPLICIT double precision (A-H,O-Z)
C
      PARAMETER (b01=-0.177 317 424 732 13 D-02,
     1           b02=-0.178 348 622 923 58 D-01,
     2           b03=-0.459 960 136 963 65 D-01,
     3           b04=-0.575 812 590 834 32 D-01,
     4           b05=-0.503 252 787 279 30 D-01,
     5           b06=-0.330 326 416 702 03 D-04*2.D0,
     6           b07=-0.189 489 875 163 15 D-03*2.D0,
     7           b08=-0.393 927 772 433 55 D-02*2.D0,
     8           b09=-0.437 972 956 505 73 D-01*2.D0,
     9           b10=-0.266 745 479 140 87 D-04*2.D0,
     1           b11= 0.204 817 376 923 09 D-07*3.D0,
     2           b12= 0.438 706 672 844 35 D-06*3.D0,
     3           b13=-0.322 776 772 385 70 D-04*3.D0,
     4           b14=-0.150 339 245 421 48 D-02*3.D0,
     5           b15=-0.406 682 535 626 49 D-01*3.D0,
     6           b16=-0.788 473 095 593 67 D-09*4.D0,
     7           b17= 0.127 907 178 522 85 D-07*4.D0,
     8           b18= 0.482 253 727 185 07 D-06*4.D0,
     9           b19= 0.229 220 763 376 61 D-05*5.D0,
     1           b20=-0.167 147 664 510 61 D-10*6.D0,
     2           b21=-0.211 714 723 213 55 D-02*6.D0,
     3           b22=-0.238 957 419 341 04 D02*6.D0,
     4           b23=-0.590 595 643 242 70 D-17*7.D0,
     5           b24=-0.126 218 088 991 01 D-05*7.D0,
     6           b25=-0.389 468 424 357 39 D-01*7.D0,
     7           b26= 0.112 562 113 604 59 D-10*8.D0,
     8           b27=-0.823 113 408 979 98 D01*8.D0,
     9           b28= 0.198 097 128 020 88 D-07*9.D0,
     1           b29= 0.104 069 652 101 74 D-18*10.D0,
     2           b30=-0.102 347 470 959 29 D-12*10.D0,
     3           b31=-0.100 181 793 795 11 D-08*10.D0,
     4           b32=-0.808 829 086 469 85 D-10*16.D0,
     5           b33= 0.106 930 318 794 09 D0*16.D0,
     6           b34=-0.336 622 505 741 71 D0*18.D0,
     7           b35= 0.891 858 453 554 21 D-24*20.D0,
     8           b36= 0.306 293 168 762 32 D-12*20.D0,
     9           b37=-0.420 024 676 982 08 D-05*20.D0,
     1           b38=-0.590 560 296 856 39 D-25*21.D0,
     2           b39= 0.378 269 476 134 57 D-05*22.D0,
     3           b40=-0.127 686 089 346 81 D-14*23.D0,
     4           b41= 0.730 876 105 950 61 D-28*24.D0,
     5           b42= 0.554 147 153 507 78 D-16*24.D0,
     6           b43=-0.943 697 072 412 10 D-06*24.D0)
      PARAMETER (a06=-0.330 326 416 702 03 D-04*2.D0,
     6           a07=-0.189 489 875 163 15 D-03*2.D0,
     7           a08=-0.393 927 772 433 55 D-02*2.D0,
     8           a09=-0.437 972 956 505 73 D-01*2.D0,
     9           a10=-0.266 745 479 140 87 D-04*2.D0,
     1           a11= 0.204 817 376 923 09 D-07*6.D0,
     2           a12= 0.438 706 672 844 35 D-06*6.D0,
     3           a13=-0.322 776 772 385 70 D-04*6.D0,
     4           a14=-0.150 339 245 421 48 D-02*6.D0,
     5           a15=-0.406 682 535 626 49 D-01*6.D0,
     6           a16=-0.788 473 095 593 67 D-09*12.D0,
     7           a17= 0.127 907 178 522 85 D-07*12.D0,
     8           a18= 0.482 253 727 185 07 D-06*12.D0,
     9           a19= 0.229 220 763 376 61 D-05*20.D0,
     1           a20=-0.167 147 664 510 61 D-10*30.D0,
     2           a21=-0.211 714 723 213 55 D-02*30.D0,
     3           a22=-0.238 957 419 341 04 D02*30.D0,
     4           a23=-0.590 595 643 242 70 D-17*42.D0,
     5           a24=-0.126 218 088 991 01 D-05*42.D0,
     6           a25=-0.389 468 424 357 39 D-01*42.D0,
     7           a26= 0.112 562 113 604 59 D-10*56.D0,
     8           a27=-0.823 113 408 979 98 D01*56.D0,
     9           a28= 0.198 097 128 020 88 D-07*72.D0,
     1           a29= 0.104 069 652 101 74 D-18*90.D0,
     2           a30=-0.102 347 470 959 29 D-12*90.D0,
     3           a31=-0.100 181 793 795 11 D-08*90.D0,
     4           a32=-0.808 829 086 469 85 D-10*240.D0,
     5           a33= 0.106 930 318 794 09 D0*240.D0,
     6           a34=-0.336 622 505 741 71 D0*306.D0,
     7           a35= 0.891 858 453 554 21 D-24*380.D0,
     8           a36= 0.306 293 168 762 32 D-12*380.D0,
     9           a37=-0.420 024 676 982 08 D-05*380.D0,
     1           a38=-0.590 560 296 856 39 D-25*420.D0,
     2           a39= 0.378 269 476 134 57 D-05*462.D0,
     3           a40=-0.127 686 089 346 81 D-14*506.D0,
     4           a41= 0.730 876 105 950 61 D-28*552.D0,
     5           a42= 0.554 147 153 507 78 D-16*552.D0,
     6           a43=-0.943 697 072 412 10 D-06*552.D0)
      PARAMETER (c03=-0.560 879 112 830 20 D-02*30.D0,
     3           c04= 0.714 527 380 814 55 D-01*20.D0,
     4           c05=-0.407 104 982 239 28 D0*12.D0,
     5           c06= 0.142 408 191 714 44 D01*6.D0,
     6           c07=-0.438 395 113 194 50 D01*2.D0,
     7           c08=-0.284 086 324 607 72 D0*2.D0,
     8           c09= 0.212 684 637 533 07 D-01*6.D0)
      PARAMETER (d03=-0.459 960 136 963 65 D-01*2.D0,
     3           d04=-0.575 812 590 834 32 D-01*6.D0,
     4           d05=-0.503 252 787 279 30 D-01*30.D0,
     6           d07=-0.189 489 875 163 15 D-03*2.D0,
     7           d08=-0.393 927 772 433 55 D-02*12.D0,
     8           d09=-0.437 972 956 505 73 D-01*42.D0,
     9           d10=-0.266 745 479 140 87 D-04*1260.D0,
     3           d13=-0.322 776 772 385 70 D-04*6.D0,
     4           d14=-0.150 339 245 421 48 D-02*30.D0,
     5           d15=-0.406 682 535 626 49 D-01*1190.D0,
     7           d17= 0.127 907 178 522 85 D-07*2.D0,
     8           d18= 0.482 253 727 185 07 D-06*6.D0,
     9           d19= 0.229 220 763 376 61 D-05*42.D0,
     1           d20=-0.167 147 664 510 61 D-10*6.D0,
     2           d21=-0.211 714 723 213 55 D-02*240.D0,
     3           d22=-0.238 957 419 341 04 D02*1190.D0,
     5           d24=-0.126 218 088 991 01 D-05*110.D0,
     6           d25=-0.389 468 424 357 39 D-01*600.D0,
     7           d26= 0.112 562 113 604 59 D-10*56.D0,
     8           d27=-0.823 113 408 979 98 D01*1260.D0,
     9           d28= 0.198 097 128 020 88 D-07*156.D0,
     1           d29= 0.104 069 652 101 74 D-18*12.D0,
     2           d30=-0.102 347 470 959 29 D-12*90.D0,
     3           d31=-0.100 181 793 795 11 D-08*182.D0,
     4           d32=-0.808 829 086 469 85 D-10*812.D0,
     5           d33= 0.106 930 318 794 09 D0*2450.D0,
     6           d34=-0.336 622 505 741 71 D0*3192.D0,
     7           d35= 0.891 858 453 554 21 D-24*380.D0,
     8           d36= 0.306 293 168 762 32 D-12*1190.D0,
     9           d37=-0.420 024 676 982 08 D-05*2256.D0,
     1           d38=-0.590 560 296 856 39 D-25*420.D0,
     2           d39= 0.378 269 476 134 57 D-05*2756.D0,
     3           d40=-0.127 686 089 346 81 D-14*1482.D0,
     4           d41= 0.730 876 105 950 61 D-28*650.D0,
     5           d42= 0.554 147 153 507 78 D-16*1560.D0,
     6           d43=-0.943 697 072 412 10 D-06*3306.D0)
      PARAMETER (e02=-0.178 348 622 923 58 D-01,
     2           e03=-0.459 960 136 963 65 D-01*2.D0,
     3           e04=-0.575 812 590 834 32 D-01*3.D0,
     4           e05=-0.503 252 787 279 30 D-01*6.D0,
     5           e06=-0.330 326 416 702 03 D-04*2.D0,
     6           e07=-0.189 489 875 163 15 D-03*4.D0,
     7           e08=-0.393 927 772 433 55 D-02*8.D0,
     8           e09=-0.437 972 956 505 73 D-01*14.D0,
     9           e10=-0.266 745 479 140 87 D-04*72.D0,
     2           e12= 0.438 706 672 844 35 D-06*3.D0,
     3           e13=-0.322 776 772 385 70 D-04*9.D0,
     4           e14=-0.150 339 245 421 48 D-02*18.D0,
     5           e15=-0.406 682 535 626 49 D-01*105.D0,
     6           e16=-0.788 473 095 593 67 D-09*4.D0,
     7           e17= 0.127 907 178 522 85 D-07*8.D0,
     8           e18= 0.482 253 727 185 07 D-06*12.D0,
     9           e19= 0.229 220 763 376 61 D-05*35.D0,
     1           e20=-0.167 147 664 510 61 D-10*18.D0,
     2           e21=-0.211 714 723 213 55 D-02*96.D0,
     3           e22=-0.238 957 419 341 04 D02*210.D0,
     5           e24=-0.126 218 088 991 01 D-05*77.D0,
     6           e25=-0.389 468 424 357 39 D-01*175.D0,
     7           e26= 0.112 562 113 604 59 D-10*64.D0,
     8           e27=-0.823 113 408 979 98 D01*288.D0,
     9           e28= 0.198 097 128 020 88 D-07*117.D0,
     1           e29= 0.104 069 652 101 74 D-18*40.D0,
     2           e30=-0.102 347 470 959 29 D-12*100.D0,
     3           e31=-0.100 181 793 795 11 D-08*140.D0,
     4           e32=-0.808 829 086 469 85 D-10*464.D0,
     5           e33= 0.106 930 318 794 09 D0*800.D0,
     6           e34=-0.336 622 505 741 71 D0*1026.D0,
     7           e35= 0.891 858 453 554 21 D-24*400.D0,
     8           e36= 0.306 293 168 762 32 D-12*700.D0,
     9           e37=-0.420 024 676 982 08 D-05*960.D0,
     1           e38=-0.590 560 296 856 39 D-25*441.D0,
     2           e39= 0.378 269 476 134 57 D-05*1166.D0,
     3           e40=-0.127 686 089 346 81 D-14*897.D0,
     4           e41= 0.730 876 105 950 61 D-28*624.D0,
     5           e42= 0.554 147 153 507 78 D-16*960.D0,
     6           e43=-0.943 697 072 412 10 D-06*1392.D0)
C
      pi= p
      tau=540.d0/t
      tauinv= 1.d0/tau
C
      taun= tau - 0.5d0
      taun2= taun*taun
      taun3= taun2*taun
      taun4= taun3*taun
      taun6= taun4*taun2
      taun7= taun6*taun
      taun9= taun7*taun2
      taun10= taun9*taun
      taun13= taun10*taun3
      taun14= taun13*taun
      taun15= taun14*taun
      taun19= taun15*taun4
      taun20= taun19*taun
      taun29= taun20*taun9
C
      pi2= pi*pi
      pi4= pi2*pi2
C
      gap2= b01 + taun*(b02 + taun*(b03 + taun*(b04 + taun3*b05))) +
     3     pi* (taun*(b06 + taun*(b07 + taun2*(b08 +
     4         taun3*(b09 + taun29*b10)))) +
     5     pi* (b11 + taun*(b12 + taun2*(b13 + taun3*(b14 +
     6         taun29*b15))) +
     7     pi* (taun*(b16 + taun*(b17 + taun*b18)) +
     8     pi* (taun7*b19 +
     9     pi* (taun3*(b20 + taun13*(b21 + taun19*b22)) +
     1     pi* (b23 + taun10*taun*(b24 + taun14*b25) +
     2     pi* (taun7*(taun*b26 + taun29*b27) +
     3     pi* (taun13*b28 +
     4     pi* (taun4*(b29 + taun6*(b30 + taun4*b31)) +
     5     pi4*pi2*taun20*(taun9*b32 + taun20*taun10*(b33 +
     6         pi2*taun7*b34) +
     7     pi4*(b35 + taun15*(b36 + taun13*b37) +
     8     pi* (taun*b38 +
     9     pi* (taun19*(taun14*b39 +
     1          pi* b40) +
     2     pi2*(taun6*(b41 + taun14*(b42 +
     3     taun14*taun4*b43))) ))) ))) ))) ))) )
C
      gapp2= (taun*(a06 + taun*(a07 + taun2*(a08 +
     4     taun3*(a09 + taun29*a10)))) +
     5     pi* (a11 + taun*(a12 + taun2*(a13 + taun3*(a14 +
     6         taun29*a15))) +
     7     pi* (taun*(a16 + taun*(a17 + taun*a18)) +
     8     pi* (taun7*a19 +
     9     pi* (taun3*(a20 + taun13*(a21 + taun19*a22)) +
     1     pi* (a23 + taun10*taun*(a24 + taun14*a25) +
     2     pi* (taun7*(taun*a26 + taun29*a27) +
     3     pi* (taun13*a28 +
     4     pi* (taun4*(a29 + taun6*(a30 + taun4*a31)) +
     5     pi4*pi2*taun20*(taun9*a32 + taun20*taun10*(a33 +
     6         pi2*taun7*a34) +
     7     pi4*(a35 + taun15*(a36 + taun13*a37) +
     8     pi* (taun*a38 +
     9     pi* (taun19*(taun14*a39 +
     1          pi* a40) +
     2     pi2*(taun6*(a41 + taun14*(a42 +
     3     taun14*taun4*a43))) ))) ))) ))) ))) )
C
      gatt2= c08 + tau*c09 +
     1     tauinv*tauinv*tauinv*(c07 + tauinv*(c06 +
     2     tauinv*(c05 + tauinv*(c04 + tauinv*c03)))) +
     3     pi* (d03 + taun*(d04 + taun3*d05) +
     3     pi* (d07 + taun2*(d08 + taun3*(d09 + taun29*d10)) +
     5     pi* (taun*(d13 + taun3*(d14 + taun29*d15)) +
     7     pi* (d17 + taun*d18 +
     8     pi* (taun3*taun2*d19 +
     9     pi* (taun*(d20 + taun13*(d21 + taun19*d22)) +
     1     pi* (taun9*(d24 + taun14*d25) +
     2     pi* (taun3*taun2*(taun*d26 + taun29*d27) +
     3     pi* (taun10*taun*d28 +
     4     pi* (taun2*(d29 + taun6*(d30 + taun4*d31)) +
     5     pi4*pi2*taun15*taun3*(taun9*d32 + taun20*taun10*(d33 +
     6         pi2*taun7*d34) +
     7     pi4*(d35 + taun15*(d36 + taun13*d37) +
     8     pi* (taun*d38 +
     9     pi* (taun19*(taun14*d39 + pi* d40) +
     2     pi2*(taun6*(d41 + taun14*(d42 +
     3     taun14*taun4*d43))) ))) ))) ))) ))) ))
C
      gapt2= e02 + taun*(e03 + taun*(e04 + taun3*e05)) +
     3     pi* (e06 + taun*(e07 + taun2*(e08 +
     4         taun3*(e09 + taun29*e10))) +
     5     pi* (e12 + taun2*(e13 + taun3*(e14 + taun29*e15)) +
     7     pi* (e16 + taun*(e17 + taun*e18) +
     8     pi* (taun6*e19 +
     9     pi* (taun2*(e20 + taun13*(e21 + taun19*e22)) +
     1     pi* (taun10*(e24 + taun14*e25) +
     2     pi* (taun6*(taun*e26 + taun29*e27) +
     3     pi* (taun6*taun6*e28 +
     4     pi* (taun3*(e29 + taun6*(e30 + taun4*e31)) +
     5     pi4*pi2*taun19*(taun9*e32 + taun20*taun10*(e33 +
     6         pi2*taun7*e34) +
     7     pi4*(e35 + taun15*(e36 + taun13*e37) +
     8     pi* (taun*e38 +
     9     pi* (taun19*(taun14*e39 + pi* e40) +
     2     pi2*(taun6*(e41 + taun14*(e42 +
     3     taun14*taun4*e43))) ))) ))) ))) ))) )
C
      wpt2n= dsqrt(0.461526d3*t*(1.d0 + 2.d0*pi*gap2 + pi2*gap2*gap2)/
     1      ((1.d0 - pi2*gapp2) + (1.d0 + pi*gap2 - tau*pi*gapt2)*
     2      (1.d0 + pi*gap2 - tau*pi*gapt2)/(tau*tau*gatt2)))
C
      END
C
c
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is pressure in MPa.
C  Second argument is temperature in K.
C  Returns specific entropy in kJ/(kg*K)
C
      FUNCTION SPT2N(P,T)
C
      IMPLICIT double precision (A-H,O-Z)
C
      parameter(a04=-0.44448764333452D+01,
     5          a05=-0.22926624714607D+02,
     6          a06=-0.43051902051180D+02,
     7          a07=-0.75253615672206D+02,
     8          a08=-0.82325284089205D-02,
     9          a09=-0.94450864454513D-01,
     &          a10=-0.39270508365637D+01,
     1          a11=-0.76407372741773D+02,
     2          a12=-0.23932578946761D+00,
     4          a14= 0.10933624938123D-03,
     5          a15=-0.24133119369638D-01,
     6          a16=-0.22480892468696D+01,
     7          a17=-0.35474272584198D+03,
     8          a18=-0.19650645031516D-06,
     9          a19= 0.63755087552931D-05,
     &          a20= 0.36056766658235D-03,
     1          a21= 0.39989127290421D-02,
     2          a22=-0.12497164867770D-07,
     3          a23=-0.84423037834823D+01,
     4          a24=-0.20843876702652D+06,
     6          a26=-0.34602240265361D-02,
     7          a27=-0.24266223542696D+03,
     8          a28= 0.22442547762780D-07,
     9          a29=-0.73850273699100D+05,
     &          a30= 0.64181736525088D-04,
     1          a31= 0.10374663655276D-15,
     2          a32=-0.25507450196258D-09,
     3          a33=-0.34954795937691D-05,
     4          a34=-0.58458099253864D-06,
     5          a35= 0.13324803024175D+04,
     6          a36=-0.47819819876448D+04,
     7          a37= 0.44454513380586D-20,
     8          a38= 0.26717467330171D-08,
     9          a39=-0.50246518510642D-01,
     &          a40=-0.30908182839692D-21,
     1          a41= 0.49965138936998D-01,
     2          a42=-0.12410752785137D-10,
     3          a43= 0.47359492924764D-24,
     4          a44= 0.55242716940682D-12,
     5          a45=-0.13641135821518D-01)
C
      parameter(b1= 0.44734647439699D+01,
     2          b3= 0.15531617605684D-01,
     3          b4=-0.16488648197891D+00,
     4          b5= 0.75155813613186D+00,
     5          b6=-0.19717524926760D+01,
     6          b7= 0.40466148602441D+01,
     7          b8=-0.13111322505090D+00,
     8          b9= 0.19631898004417D-01,
     8          b0=-0.461526D+00)
C
      parameter(c03=-0.81836601766922D-03,
     4          c04=-0.82312526543430D-02,
     5          c05=-0.21228356217229D-01,
     6          c06=-0.26575248179740D-01,
     7          c07=-0.23226424590187D-01,
     8          c08=-0.15245422979482D-04,
     9          c09=-0.87454504124548D-04,
     &          c10=-0.18180790910017D-02,
     1          c11=-0.20213590672427D-01,
     2          c12=-0.12310997400597D-04,
     3          c13= 0.94528544701805D-08,
     4          c14= 0.20247453589116D-06,
     5          c15=-0.14896987265208D-04,
     6          c16=-0.69385470582394D-03,
     7          c17=-0.18769456393755D-01,
     8          c18=-0.36390083391697D-09,
     9          c19= 0.59032488474936D-08,
     &          c20= 0.22257263369282D-06,
     1          c21= 0.10579134203815D-05,
     2          c22=-0.77142993010924D-11,
     3          c23=-0.97711849345858D-03,
     4          c24=-0.11028506191879D+02,
     5          c25=-0.27257524484323D-17,
     6          c26=-0.58252929739665D-06,
     7          c27=-0.17974980401997D-01,
     8          c28= 0.51950342043472D-11,
     9          c29=-0.37988823919289D+01,
     &          c30= 0.91426975106963D-08,
     1          c31= 0.48030850255907D-19,
     2          c32=-0.47236018881958D-13,
     3          c33=-0.46236502563082D-09,
     4          c34=-0.37329565296209D-10,
     5          c35= 0.49351122311761D-01,
     6          c36=-0.15536003858495D+00,
     7          c37= 0.41161586463506D-24,
     8          c38= 0.14136226100620D-12,
     9          c39=-0.19385230906884D-05,
     &          c40=-0.27255893156695D-25,
     1          c41= 0.17458119824248D-05,
     2          c42=-0.58930450071876D-15,
     3          c43= 0.33731832567496D-28,
     4          c44= 0.25575331916983D-16,
     5          c45=-0.43554073504207D-06)
C
      tau = 540.d0/t
      ti= t*0.18518518518518D-02
C
      targ = tau - 0.5d0
      targ2 = targ*targ
      targ3 = targ2*targ
      targ4 = targ3*targ
      targ6 = targ4*targ2
      targ7 = targ6*targ
      targ9 = targ7*targ2
      targ10 = targ9*targ
      targ13 = targ10*targ3
      targ14 = targ13*targ
      targ15 = targ14*targ
      targ18 = targ15*targ3
      targ19 = targ18*targ
      targ20 = targ19*targ
      targ28 = targ19*targ9
      targ29 = targ28*targ
      targ30 = targ29*targ
C
      pi = p
      pi2 = pi*pi
      pi4 = pi2*pi2
      pi6 = pi4*pi2
C
      p2t6 = pi2*targ6
      p2t7 = p2t6*targ
C
      s1=((((((((((((((((targ18*a45 + a44)*targ14 + a43)*
     &   p2t6 + (pi*a42 + targ14*a41)*targ19)*pi +
     &   targ*a40)*pi + (targ13*a39 + a38)*targ15 + a37)*
     &   pi4 + (p2t7*a36 + a35)*targ30 + targ9*a34)*
     &   pi6*targ19 + targ13*a33 + targ9*a32 +
     &   targ3*a31)*pi + targ10*targ2*a30)*pi + (targ28*a29 +
     &   a28)*targ7)*pi + (targ14*a27 + a26)*targ10)*pi + (targ19*
     &   a24 + a23)*targ15 + targ2*a22)*pi + targ6*
     &   a21)*pi + targ2*a20 + targ*a19 + a18)*pi + ((targ29*
     &   a17 + a16)*targ3 + a15)*targ2 + a14)*pi + (((targ29*
     &   a12 + a11)*targ3 + a10)*targ2 + a09)*targ + a08)*pi + (
     &   targ3*a07 + a06)*targ2 + targ*a05 + a04)*pi/t
     &   +((((b3*ti + b4)*ti + b5)*ti + b6)*ti + b7)*ti
     &   + (b9*tau + b8)*tau*tau + b0*dlog(pi) + b1
C
      s2=((((((((((((((((targ18*c45 + c44)*targ14 + c43)*
     &   p2t6 + (pi*c42 + targ14*c41)*targ19)*pi + targ*c40)*
     &   pi + (targ13*c39 + c38)*targ15 + c37)*pi4 + (p2t7*
     &   c36 + c35)*targ30 + targ9*c34)*pi6*targ20 +
     &   targ14*c33 + targ10*c32 + targ4*c31)*pi + targ13*c30)*pi +
     &   (targ29*c29 + targ*c28)*targ7)*pi + (targ14*c27 +
     &   c26)*targ10*targ + c25)*pi + (targ19*c24 + c23)*targ14*targ2 +
     &   targ3*c22)*pi + targ7*c21)*pi + targ3*c20 + targ2*c19 + targ*
     &   c18)*pi + ((targ29*c17 + c16)*targ3 + c15)*targ3 + targ*c14 +
     &   c13)*pi + (((targ29*c12 + c11)*targ3 + c10)*targ2 + c09)*
     &   targ2 + targ*c08)*pi + (targ3*c07 + c06)*targ3 + targ2*c05 +
     &   targ*c04 + c03)*pi
C
      spt2n=s1-s2
C
      END
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is pressure in MPa.
C  Second argument is temperature in K.
C  Returns specific isobaric heat capacity in kJ/(kg*K)
C
      FUNCTION CPPT2N(P,T)
C
      IMPLICIT double precision (A-H,O-Z)
C
      PARAMETER (a03=-0.560 879 112 830 20 D-02*30.D0,
     3           a04= 0.714 527 380 814 55 D-01*20.D0,
     4           a05=-0.407 104 982 239 28 D0*12.D0,
     5           a06= 0.142 408 191 714 44 D01*6.D0,
     6           a07=-0.438 395 113 194 50 D01*2.D0,
     7           a08=-0.284 086 324 607 72 D0*2.D0,
     8           a09= 0.212 684 637 533 07 D-01*6.D0)
      PARAMETER (b03=-0.459 960 136 963 65 D-01*2.D0,
     3           b04=-0.575 812 590 834 32 D-01*6.D0,
     4           b05=-0.503 252 787 279 30 D-01*30.D0,
     6           b07=-0.189 489 875 163 15 D-03*2.D0,
     7           b08=-0.393 927 772 433 55 D-02*12.D0,
     8           b09=-0.437 972 956 505 73 D-01*42.D0,
     9           b10=-0.266 745 479 140 87 D-04*1260.D0,
     3           b13=-0.322 776 772 385 70 D-04*6.D0,
     4           b14=-0.150 339 245 421 48 D-02*30.D0,
     5           b15=-0.406 682 535 626 49 D-01*1190.D0,
     7           b17= 0.127 907 178 522 85 D-07*2.D0,
     8           b18= 0.482 253 727 185 07 D-06*6.D0,
     9           b19= 0.229 220 763 376 61 D-05*42.D0,
     1           b20=-0.167 147 664 510 61 D-10*6.D0,
     2           b21=-0.211 714 723 213 55 D-02*240.D0,
     3           b22=-0.238 957 419 341 04 D02*1190.D0,
     5           b24=-0.126 218 088 991 01 D-05*110.D0,
     6           b25=-0.389 468 424 357 39 D-01*600.D0,
     7           b26= 0.112 562 113 604 59 D-10*56.D0,
     8           b27=-0.823 113 408 979 98 D01*1260.D0,
     9           b28= 0.198 097 128 020 88 D-07*156.D0,
     1           b29= 0.104 069 652 101 74 D-18*12.D0,
     2           b30=-0.102 347 470 959 29 D-12*90.D0,
     3           b31=-0.100 181 793 795 11 D-08*182.D0,
     4           b32=-0.808 829 086 469 85 D-10*812.D0,
     5           b33= 0.106 930 318 794 09 D0*2450.D0,
     6           b34=-0.336 622 505 741 71 D0*3192.D0,
     7           b35= 0.891 858 453 554 21 D-24*380.D0,
     8           b36= 0.306 293 168 762 32 D-12*1190.D0,
     9           b37=-0.420 024 676 982 08 D-05*2256.D0,
     1           b38=-0.590 560 296 856 39 D-25*420.D0,
     2           b39= 0.378 269 476 134 57 D-05*2756.D0,
     3           b40=-0.127 686 089 346 81 D-14*1482.D0,
     4           b41= 0.730 876 105 950 61 D-28*650.D0,
     5           b42= 0.554 147 153 507 78 D-16*1560.D0,
     6           b43=-0.943 697 072 412 10 D-06*3306.D0)
C
      PARAMETER ( tn= 540.0D0,
     1            tnq= 1.D0/tn,
     2            r= 0.461526D0)
C
      tau=tn/t
      tauinv= t*tnq
C
      taun= tau - 0.5d0
      taun2= taun*taun
      taun3= taun2*taun
      taun4= taun3*taun
      taun6= taun4*taun2
      taun7= taun6*taun
      taun9= taun7*taun2
      taun10= taun9*taun
      taun13= taun10*taun3
      taun14= taun13*taun
      taun15= taun14*taun
      taun19= taun15*taun4
      taun20= taun19*taun
      taun29= taun20*taun9
C
      pi= p
      pi2= pi*pi
      pi4= pi2*pi2
C
      cppt2n= (a08 + tau*a09 +
     1     tauinv*tauinv*tauinv*(a07 + tauinv*(a06 +
     2     tauinv*(a05 + tauinv*(a04 + tauinv*a03)))) +
     3     pi* (b03 + taun*(b04 + taun3*b05) +
     3     pi* (b07 + taun2*(b08 + taun3*(b09 + taun29*b10)) +
     5     pi* (taun*(b13 + taun3*(b14 + taun29*b15)) +
     7     pi* (b17 + taun*b18 +
     8     pi* (taun3*taun2*b19 +
     9     pi* (taun*(b20 + taun13*(b21 + taun19*b22)) +
     1     pi* (taun9*(b24 + taun14*b25) +
     2     pi* (taun3*taun2*(taun*b26 + taun29*b27) +
     3     pi* (taun10*taun*b28 +
     4     pi* (taun2*(b29 + taun6*(b30 + taun4*b31)) +
     5     pi4*pi2*taun15*taun3*(taun9*b32 + taun20*taun10*(b33 +
     6         pi2*taun7*b34) +
     7     pi4*(b35 + taun15*(b36 + taun13*b37) +
     8     pi* (taun*b38 +
     9     pi* (taun19*(taun14*b39 + pi* b40) +
     2     pi2*(taun6*(b41 + taun14*(b42 +
     3     taun14*taun4*b43))) ))) ))) ))) ))) )))*tau*tau*(-r)
C
      END
C
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is pressure in MPa.
C  Second argument is temperature in K.
C  Returns specific enthalpy in kJ/kg
C
      FUNCTION HPT2N(P,T)
C
      IMPLICIT double precision (A-H,O-Z)
C
      parameter(a04=-0.44448764333452D+01,
     5          a05=-0.22926624714607D+02,
     6          a06=-0.43051902051180D+02,
     7          a07=-0.75253615672206D+02,
     8          a08=-0.82325284089205D-02,
     9          a09=-0.94450864454513D-01,
     &          a10=-0.39270508365637D+01,
     1          a11=-0.76407372741773D+02,
     2          a12=-0.23932578946761D+00,
     4          a14= 0.10933624938123D-03,
     5          a15=-0.24133119369638D-01,
     6          a16=-0.22480892468696D+01,
     7          a17=-0.35474272584198D+03,
     8          a18=-0.19650645031516D-06,
     9          a19= 0.63755087552931D-05,
     &          a20= 0.36056766658235D-03,
     1          a21= 0.39989127290421D-02,
     2          a22=-0.12497164867770D-07,
     3          a23=-0.84423037834823D+01,
     4          a24=-0.20843876702652D+06,
     6          a26=-0.34602240265361D-02,
     7          a27=-0.24266223542696D+03,
     8          a28= 0.22442547762780D-07,
     9          a29=-0.73850273699100D+05,
     &          a30= 0.64181736525088D-04,
     1          a31= 0.10374663655276D-15,
     2          a32=-0.25507450196258D-09,
     3          a33=-0.34954795937691D-05,
     4          a34=-0.58458099253864D-06,
     5          a35= 0.13324803024175D+04,
     6          a36=-0.47819819876448D+04,
     7          a37= 0.44454513380586D-20,
     8          a38= 0.26717467330171D-08,
     9          a39=-0.50246518510642D-01,
     &          a40=-0.30908182839692D-21,
     1          a41= 0.49965138936998D-01,
     2          a42=-0.12410752785137D-10,
     3          a43= 0.47359492924764D-24,
     4          a44= 0.55242716940682D-12,
     5          a45=-0.13641135821518D-01)
C
      parameter(b2= 0.25138371504395D+04,
     2          b3= 0.69892279225579D+01,
     3          b4=-0.71230960214888D+02,
     4          b5= 0.30438104513340D+03,
     5          b6=-0.70983089736335D+03,
     6          b7= 0.10925860122659D+04,
     7          b8=-0.14160228305498D+03,
     8          b9= 0.15901837383578D+02)
C
      tau= 540.d0/t
      ti = t*0.18518518518518D-02
C
      targ = tau - 0.5d0
      targ2 = targ*targ
      targ3 = targ2*targ
      targ4 = targ3*targ
      targ6 = targ4*targ2
      targ7 = targ6*targ
      targ9 = targ7*targ2
      targ10 = targ9*targ
      targ13 = targ10*targ3
      targ14 = targ13*targ
      targ15 = targ14*targ
      targ19 = targ15*targ4
      targ20 = targ19*targ
      targ29 = targ20*targ9
C
      pi = p
      pi2 = pi*pi
      pi4 = pi2*pi2
C
      hpt2n=((((((((((((((((targ14*targ4*a45 + a44)*targ14 + a43)*
     &   pi2*targ6 + (pi*a42 + targ14*a41)*targ19)*pi +
     &   targ*a40)*pi + (targ13*a39 + a38)*targ15 + a37)*
     &   pi4 + (targ7*pi2*a36 + a35)*targ20*targ10 + targ9*a34)*
     &   pi4*pi2*targ19 + targ13*a33 + targ9*a32 +
     &   targ3*a31)*pi + targ10*targ2*a30)*pi + (targ14*targ14*a29 +
     &   a28)*targ7)*pi + (targ14*a27 + a26)*targ10)*pi + (targ19*
     &   a24 + a23)*targ15 + targ2*a22)*pi + targ6*
     &   a21)*pi + targ2*a20 + targ*a19 + a18)*pi + ((targ29*
     &   a17 + a16)*targ3 + a15)*targ2 + a14)*pi + (((targ29*
     &   a12 + a11)*targ3 + a10)*targ2 + a09)*targ + a08)*pi + (
     &   targ3*a07 + a06)*targ2 + targ*a05 + a04)*pi + ((((b3*ti +
     &   b4)*ti + b5)*ti + b6)*ti + b7)*ti*ti + b2 + (b9*tau + b8)*tau
c
      END
c
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is pressure in MPa.
C  Second argument is temperature in K.
C  Returns specific volume in m^3/kg
C
      FUNCTION VPT2N(P,T)
C
      IMPLICIT double precision (A-H,O-Z)
C
      parameter(a03=-0.81836601766921D-06,
     4          a04=-0.82312526543429D-05,
     5          a05=-0.21228356217229D-04,
     6          a06=-0.26575248179740D-04,
     7          a07=-0.23226424590187D-04,
     8          a08=-0.30490845958965D-07,
     9          a09=-0.17490900824910D-06,
     &          a10=-0.36361581820033D-05,
     1          a11=-0.40427181344853D-04,
     2          a12=-0.24621994801194D-07,
     3          a13= 0.28358563410542D-10,
     4          a14= 0.60742360767348D-09,
     5          a15=-0.44690961795625D-07,
     6          a16=-0.20815641174718D-05,
     7          a17=-0.56308369181266D-04,
     8          a18=-0.14556033356679D-11,
     9          a19= 0.23612995389975D-10,
     &          a20= 0.89029053477126D-09,
     1          a21= 0.52895671019076D-08,
     2          a22=-0.46285795806554D-13,
     3          a23=-0.58627109607515D-05,
     4          a24=-0.66171037151276D-01,
     5          a25=-0.19080267139026D-19,
     6          a26=-0.40777050817765D-08,
     7          a27=-0.12582486281398D-03,
     8          a28= 0.41560273634777D-13,
     9          a29=-0.30391059135431D-01,
     &          a30= 0.82284277596267D-10,
     1          a31= 0.48030850255907D-21,
     2          a32=-0.47236018881958D-15,
     3          a33=-0.46236502563083D-11,
     4          a34=-0.59727304473934D-12,
     5          a35= 0.78961795698818D-03,
     6          a36=-0.27964806945291D-02,
     7          a37= 0.82323172927012D-26,
     8          a38= 0.28272452201240D-14,
     9          a39=-0.38770461813767D-07,
     &          a40=-0.57237375629058D-27,
     1          a41= 0.38407863613346D-07,
     2          a42=-0.13554003516532D-16,
     3          a43= 0.80956398161990D-30,
     4          a44= 0.61380796600759D-18,
     5          a45=-0.10452977641010D-07,
     6          r=0.461526d-03)
C
      targ = 540.0d0/t - 0.5d0
      targ2 = targ*targ
      targ3 = targ2*targ
      targ4 = targ3*targ
      targ7 = targ4*targ3
      targ9 = targ7*targ2
      targ10 = targ9*targ
      targ13 = targ10*targ3
      targ14 = targ13*targ
      targ19 = targ10*targ9
      targ20 = targ19*targ
      targ29 = targ20*targ9
C
      pi = p
      pi2 = pi*pi
      pi4 = pi2*pi2
C
      vpt2n=((((((((((((((((targ14*targ4*a45 + a44)*targ14 +
     &   a43)*pi2*targ4*targ2 + (pi*a42 + targ14*a41)*targ19)*
     &   pi + targ*a40)*pi + (targ13*a39 + a38)*targ14*targ +
     &   a37)*pi4 + (targ7*pi2*a36 + a35)*targ20*targ10 + targ9*
     &   a34)*pi4*pi2*targ20 + targ14*a33 + targ10*a32 + targ4*
     &   a31)*pi + targ13*a30)*pi + (targ29*a29 + targ*a28)*
     &   targ7)*pi + (targ14*a27 + a26)*targ10*targ + a25)*pi +
     &   (targ19*a24 + a23)*targ14*targ2 + targ3*a22)*pi + targ7*
     &   a21)*pi + targ3*a20 + targ2*a19 + targ*a18)*pi + ((targ29*
     &   a17 + a16)*targ3 + a15)*targ3 + targ*a14 + a13)*pi +
     &   (((a12*targ29 + a11)*targ3 + a10)*targ2 + a09)*targ2 +
     &   targ*a08)*pi + (targ3*a07 + a06)*targ3 + targ2*a05 +
     &   targ*a04 + r/pi + a03)*t
c
      END
c
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is pressure in MPa.
C  Second argument is temperature in K.
C  Returns specific Gibbs free energy in kJ/kg.
C
      FUNCTION GPT2N(P,T)
C
      IMPLICIT double precision (A-H,O-Z)
C
      PARAMETER (a01=-0.969 276 865 002 17 D01,
     1           a02= 0.100 866 559 680 18 D02,
     2           a03=-0.560 879 112 830 20 D-02,
     3           a04= 0.714 527 380 814 55 D-01,
     4           a05=-0.407 104 982 239 28 D0,
     5           a06= 0.142 408 191 714 44 D01,
     6           a07=-0.438 395 113 194 50 D01,
     7           a08=-0.284 086 324 607 72 D0,
     8           a09= 0.212 684 637 533 07 D-01)
      PARAMETER (b01=-0.177 317 424 732 13 D-02,
     1           b02=-0.178 348 622 923 58 D-01,
     2           b03=-0.459 960 136 963 65 D-01,
     3           b04=-0.575 812 590 834 32 D-01,
     4           b05=-0.503 252 787 279 30 D-01,
     5           b06=-0.330 326 416 702 03 D-04,
     6           b07=-0.189 489 875 163 15 D-03,
     7           b08=-0.393 927 772 433 55 D-02,
     8           b09=-0.437 972 956 505 73 D-01,
     9           b10=-0.266 745 479 140 87 D-04,
     1           b11= 0.204 817 376 923 09 D-07,
     2           b12= 0.438 706 672 844 35 D-06,
     3           b13=-0.322 776 772 385 70 D-04,
     4           b14=-0.150 339 245 421 48 D-02,
     5           b15=-0.406 682 535 626 49 D-01,
     6           b16=-0.788 473 095 593 67 D-09,
     7           b17= 0.127 907 178 522 85 D-07,
     8           b18= 0.482 253 727 185 07 D-06,
     9           b19= 0.229 220 763 376 61 D-05,
     1           b20=-0.167 147 664 510 61 D-10,
     2           b21=-0.211 714 723 213 55 D-02,
     3           b22=-0.238 957 419 341 04 D02,
     4           b23=-0.590 595 643 242 70 D-17,
     5           b24=-0.126 218 088 991 01 D-05,
     6           b25=-0.389 468 424 357 39 D-01,
     7           b26= 0.112 562 113 604 59 D-10,
     8           b27=-0.823 113 408 979 98 D01,
     9           b28= 0.198 097 128 020 88 D-07,
     1           b29= 0.104 069 652 101 74 D-18,
     2           b30=-0.102 347 470 959 29 D-12,
     3           b31=-0.100 181 793 795 11 D-08,
     4           b32=-0.808 829 086 469 85 D-10,
     5           b33= 0.106 930 318 794 09 D0,
     6           b34=-0.336 622 505 741 71 D0,
     7           b35= 0.891 858 453 554 21 D-24,
     8           b36= 0.306 293 168 762 32 D-12,
     9           b37=-0.420 024 676 982 08 D-05,
     1           b38=-0.590 560 296 856 39 D-25,
     2           b39= 0.378 269 476 134 57 D-05,
     3           b40=-0.127 686 089 346 81 D-14,
     4           b41= 0.730 876 105 950 61 D-28,
     5           b42= 0.554 147 153 507 78 D-16,
     6           b43=-0.943 697 072 412 10 D-06)
C
      PARAMETER ( tn= 540.0D0,
     1            tnq= 1.D0/tn,
     2            r= 0.461526D0)
C
      tau= tn/t
      tauinv= t*tnq
C
      taun= tau - 0.5d0
      taun2= taun*taun
      taun3= taun2*taun  
      taun4= taun3*taun
      taun6= taun4*taun2
      taun7= taun6*taun
      taun9= taun7*taun2
      taun10= taun9*taun
      taun13= taun10*taun3
      taun14= taun13*taun
      taun15= taun14*taun
      taun19= taun15*taun4
      taun20= taun19*taun
      taun29= taun20*taun9
C
      pi= p
      pi2= pi*pi
      pi4= pi2*pi2
C
      gpt2n= (DLOG(pi) + a01 + tau*(a02 + tau*(a08 + tau*a09)) +
     1     tauinv*(a07 + tauinv*(a06 + tauinv*(a05 +
     2     tauinv*(a04 + tauinv*a03)))) +
     3     pi* (b01 + taun*(b02 + taun*(b03 + taun*(b04 +
     2         taun3*b05))) +
     3     pi* (taun*(b06 + taun*(b07 + taun2*(b08 +
     4         taun3*(b09 + taun29*b10)))) +
     5     pi* (b11 + taun*(b12 + taun2*(b13 + taun3*(b14 +
     6         taun29*b15))) +
     7     pi* (taun*(b16 + taun*(b17 + taun*b18)) +
     8     pi* (taun7*b19 +
     9     pi* (taun3*(b20 + taun13*(b21 + taun19*b22)) +
     1     pi* (b23 + taun10*taun*(b24 + taun14*b25) +
     2     pi* (taun7*(taun*b26 + taun29*b27) +
     3     pi* (taun13*b28 +
     4     pi* (taun4*(b29 + taun6*(b30 + taun4*b31)) +
     5     pi4*pi2*taun20*(taun9*b32 + taun20*taun10*(b33 +
     6         pi2*taun7*b34) +
     7     pi4*(b35 + taun15*(b36 + taun13*b37) +
     8     pi* (taun*b38 +
     9     pi* (taun19*(taun14*b39 +
     1          pi* b40) +
     2     pi2*(taun6*(b41 + taun14*(b42 +
     3     taun14*taun4*b43))) ))) ))) ))) ))) ))) * t*r
C
      END
C
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is pressure in MPa.
C  Second argument is temperature in K.
C  Returns specific internal energy in kJ/kg
C
      FUNCTION UPT2N(P,T)
C
      IMPLICIT double precision (A-H,O-Z)
C
      PARAMETER (a02= 0.100 866 559 680 18 D02,
     2           a03= 0.560 879 112 830 20 D-02*5.D0,
     3           a04=-0.714 527 380 814 55 D-01*4.D0,
     4           a05= 0.407 104 982 239 28 D0*3.D0,
     5           a06=-0.142 408 191 714 44 D01*2.D0,
     6           a07= 0.438 395 113 194 50 D01,
     7           a08=-0.284 086 324 607 72 D0*2.D0,
     8           a09= 0.212 684 637 533 07 D-01*3.D0)
      PARAMETER (b02=-0.178 348 622 923 58 D-01,
     2           b03=-0.459 960 136 963 65 D-01*2.D0,
     3           b04=-0.575 812 590 834 32 D-01*3.D0,
     4           b05=-0.503 252 787 279 30 D-01*6.D0,
     5           b06=-0.330 326 416 702 03 D-04,
     6           b07=-0.189 489 875 163 15 D-03*2.D0,
     7           b08=-0.393 927 772 433 55 D-02*4.D0,
     8           b09=-0.437 972 956 505 73 D-01*7.D0,
     9           b10=-0.266 745 479 140 87 D-04*36.D0,
     2           b12= 0.438 706 672 844 35 D-06,
     3           b13=-0.322 776 772 385 70 D-04*3.D0,
     4           b14=-0.150 339 245 421 48 D-02*6.D0,
     5           b15=-0.406 682 535 626 49 D-01*35.D0,
     6           b16=-0.788 473 095 593 67 D-09,
     7           b17= 0.127 907 178 522 85 D-07*2.D0,
     8           b18= 0.482 253 727 185 07 D-06*3.D0,
     9           b19= 0.229 220 763 376 61 D-05*7.D0,
     1           b20=-0.167 147 664 510 61 D-10*3.D0,
     2           b21=-0.211 714 723 213 55 D-02*16.D0,
     3           b22=-0.238 957 419 341 04 D02*35.D0,
     5           b24=-0.126 218 088 991 01 D-05*11.D0,
     6           b25=-0.389 468 424 357 39 D-01*25.D0,
     7           b26= 0.112 562 113 604 59 D-10*8.D0,
     8           b27=-0.823 113 408 979 98 D01*36.D0,
     9           b28= 0.198 097 128 020 88 D-07*13.D0,
     1           b29= 0.104 069 652 101 74 D-18*4.D0,
     2           b30=-0.102 347 470 959 29 D-12*10.D0,
     3           b31=-0.100 181 793 795 11 D-08*14.D0,
     4           b32=-0.808 829 086 469 85 D-10*29.D0,
     5           b33= 0.106 930 318 794 09 D0*50.D0,
     6           b34=-0.336 622 505 741 71 D0*57.D0,
     7           b35= 0.891 858 453 554 21 D-24*20.D0,
     8           b36= 0.306 293 168 762 32 D-12*35.D0,
     9           b37=-0.420 024 676 982 08 D-05*48.D0,
     1           b38=-0.590 560 296 856 39 D-25*21.D0,
     2           b39= 0.378 269 476 134 57 D-05*53.D0,
     3           b40=-0.127 686 089 346 81 D-14*39.D0,
     4           b41= 0.730 876 105 950 61 D-28*26.D0,
     5           b42= 0.554 147 153 507 78 D-16*40.D0,
     6           b43=-0.943 697 072 412 10 D-06*58.D0)
      PARAMETER (c01=-0.177 317 424 732 13 D-02,
     1           c02=-0.178 348 622 923 58 D-01,
     2           c03=-0.459 960 136 963 65 D-01,
     3           c04=-0.575 812 590 834 32 D-01,
     4           c05=-0.503 252 787 279 30 D-01,
     5           c06=-0.330 326 416 702 03 D-04*2.D0,
     6           c07=-0.189 489 875 163 15 D-03*2.D0,
     7           c08=-0.393 927 772 433 55 D-02*2.D0,
     8           c09=-0.437 972 956 505 73 D-01*2.D0,
     9           c10=-0.266 745 479 140 87 D-04*2.D0,
     1           c11= 0.204 817 376 923 09 D-07*3.D0,
     2           c12= 0.438 706 672 844 35 D-06*3.D0,
     3           c13=-0.322 776 772 385 70 D-04*3.D0,
     4           c14=-0.150 339 245 421 48 D-02*3.D0,
     5           c15=-0.406 682 535 626 49 D-01*3.D0,
     6           c16=-0.788 473 095 593 67 D-09*4.D0,
     7           c17= 0.127 907 178 522 85 D-07*4.D0,
     8           c18= 0.482 253 727 185 07 D-06*4.D0,
     9           c19= 0.229 220 763 376 61 D-05*5.D0,
     1           c20=-0.167 147 664 510 61 D-10*6.D0,
     2           c21=-0.211 714 723 213 55 D-02*6.D0,
     3           c22=-0.238 957 419 341 04 D02*6.D0,
     4           c23=-0.590 595 643 242 70 D-17*7.D0,
     5           c24=-0.126 218 088 991 01 D-05*7.D0,
     6           c25=-0.389 468 424 357 39 D-01*7.D0,
     7           c26= 0.112 562 113 604 59 D-10*8.D0,
     8           c27=-0.823 113 408 979 98 D01*8.D0,
     9           c28= 0.198 097 128 020 88 D-07*9.D0,
     1           c29= 0.104 069 652 101 74 D-18*10.D0,
     2           c30=-0.102 347 470 959 29 D-12*10.D0,
     3           c31=-0.100 181 793 795 11 D-08*10.D0,
     4           c32=-0.808 829 086 469 85 D-10*16.D0,
     5           c33= 0.106 930 318 794 09 D0*16.D0,
     6           c34=-0.336 622 505 741 71 D0*18.D0,
     7           c35= 0.891 858 453 554 21 D-24*20.D0,
     8           c36= 0.306 293 168 762 32 D-12*20.D0,
     9           c37=-0.420 024 676 982 08 D-05*20.D0,
     1           c38=-0.590 560 296 856 39 D-25*21.D0,
     2           c39= 0.378 269 476 134 57 D-05*22.D0,
     3           c40=-0.127 686 089 346 81 D-14*23.D0,
     4           c41= 0.730 876 105 950 61 D-28*24.D0,
     5           c42= 0.554 147 153 507 78 D-16*24.D0,
     6           c43=-0.943 697 072 412 10 D-06*24.D0)
C
      PARAMETER ( tn= 540.0D0,
     1            tnq= 1.D0/tn,
     2            r= 0.461526D0)
C
      tau= tn/t
      tauinv= t*tnq
C
      taun= tau - 0.5d0
      taun2= taun*taun
      taun3= taun2*taun
      taun4= taun3*taun
      taun6= taun4*taun2
      taun7= taun6*taun
      taun9= taun7*taun2
      taun10= taun9*taun
      taun13= taun10*taun3
      taun14= taun13*taun
      taun15= taun14*taun
      taun19= taun15*taun4
      taun20= taun19*taun
      taun29= taun20*taun9
C
      pi= p
      pi2= Pi*pi
      pi4= pi2*pi2
C
      gat2= a02 + tau*(a08 + tau*a09) +
     1     tauinv*tauinv*(a07 + tauinv*(a06 + tauinv*(a05 +
     2     tauinv*(a04 + tauinv*a03)))) +
     3     pi* (b02 + taun*(b03 + taun*(b04 + taun3*b05)) +
     3     pi* (b06 + taun*(b07 + taun2*(b08 +
     4         taun3*(b09 + taun29*b10))) +
     5     pi* (b12 + taun2*(b13 + taun3*(b14 + taun29*b15)) +
     7     pi* (b16 + taun*(b17 + taun*b18) +
     8     pi* (taun6*b19 +
     9     pi* (taun2*(b20 + taun13*(b21 + taun19*b22)) +
     1     pi* (taun10*(b24 + taun14*b25) +
     2     pi* (taun6*(taun*b26 + taun29*b27) +
     3     pi* (taun6*taun6*b28 +
     4     pi* (taun3*(b29 + taun6*(b30 + taun4*b31)) +
     5     pi4*pi2*taun19*(taun9*b32 + taun20*taun10*(b33 +
     6         pi2*taun7*b34) +
     7     pi4*(b35 + taun15*(b36 + taun13*b37) +
     8     pi* (taun*b38 +
     9     pi* (taun19*(taun14*b39 + pi* b40) +
     2     pi2*(taun6*(b41 + taun14*(b42 +
     3     taun14*taun4*b43))) ))) ))) ))) ))) ))
C
      gap2= c01 + taun*(c02 + taun*(c03 + taun*(c04 + taun3*c05))) +
     3     pi* (taun*(c06 + taun*(c07 + taun2*(c08 +
     4         taun3*(c09 + taun29*c10)))) +
     5     pi* (c11 + taun*(c12 + taun2*(c13 + taun3*(c14 +
     6         taun29*c15))) +
     7     pi* (taun*(c16 + taun*(c17 + taun*c18)) +
     8     pi* (taun7*c19 +
     9     pi* (taun3*(c20 + taun13*(c21 + taun19*c22)) +
     1     pi* (c23 + taun10*taun*(c24 + taun14*c25) +
     2     pi* (taun7*(taun*c26 + taun29*c27) +
     3     pi* (taun13*c28 +
     4     pi* (taun4*(c29 + taun6*(c30 + taun4*c31)) +
     5     pi4*pi2*taun20*(taun9*c32 + taun20*taun10*(c33 +
     6         pi2*taun7*c34) +
     7     pi4*(c35 + taun15*(c36 + taun13*c37) +
     8     pi* (taun*c38 +
     9     pi* (taun19*(taun14*c39 +
     1          pi* c40) +
     2     pi2*(taun6*(c41 + taun14*(c42 +
     3     taun14*taun4*c43))) ))) ))) ))) ))) )
C
      upt2n= r*t*(tau*gat2 - 1.D0 - pi*gap2)
C
      END
C
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is pressure in MPa.
C  Second argument is temperature in K.
C  Returns isochoric heat capacity in kJ/(kg * K)
C
      FUNCTION CVPT2N(P,T)
C
      IMPLICIT double precision (A-H,O-Z)
C
      PARAMETER (b01=-0.177 317 424 732 13 D-02,
     1           b02=-0.178 348 622 923 58 D-01,
     2           b03=-0.459 960 136 963 65 D-01,
     3           b04=-0.575 812 590 834 32 D-01,
     4           b05=-0.503 252 787 279 30 D-01,
     5           b06=-0.330 326 416 702 03 D-04*2.D0,
     6           b07=-0.189 489 875 163 15 D-03*2.D0,
     7           b08=-0.393 927 772 433 55 D-02*2.D0,
     8           b09=-0.437 972 956 505 73 D-01*2.D0,
     9           b10=-0.266 745 479 140 87 D-04*2.D0,
     1           b11= 0.204 817 376 923 09 D-07*3.D0,
     2           b12= 0.438 706 672 844 35 D-06*3.D0,
     3           b13=-0.322 776 772 385 70 D-04*3.D0,
     4           b14=-0.150 339 245 421 48 D-02*3.D0,
     5           b15=-0.406 682 535 626 49 D-01*3.D0,
     6           b16=-0.788 473 095 593 67 D-09*4.D0,
     7           b17= 0.127 907 178 522 85 D-07*4.D0,
     8           b18= 0.482 253 727 185 07 D-06*4.D0,
     9           b19= 0.229 220 763 376 61 D-05*5.D0,
     1           b20=-0.167 147 664 510 61 D-10*6.D0,
     2           b21=-0.211 714 723 213 55 D-02*6.D0,
     3           b22=-0.238 957 419 341 04 D02*6.D0,
     4           b23=-0.590 595 643 242 70 D-17*7.D0,
     5           b24=-0.126 218 088 991 01 D-05*7.D0,
     6           b25=-0.389 468 424 357 39 D-01*7.D0,
     7           b26= 0.112 562 113 604 59 D-10*8.D0,
     8           b27=-0.823 113 408 979 98 D01*8.D0,
     9           b28= 0.198 097 128 020 88 D-07*9.D0,
     1           b29= 0.104 069 652 101 74 D-18*10.D0,
     2           b30=-0.102 347 470 959 29 D-12*10.D0,
     3           b31=-0.100 181 793 795 11 D-08*10.D0,
     4           b32=-0.808 829 086 469 85 D-10*16.D0,
     5           b33= 0.106 930 318 794 09 D0*16.D0,
     6           b34=-0.336 622 505 741 71 D0*18.D0,
     7           b35= 0.891 858 453 554 21 D-24*20.D0,
     8           b36= 0.306 293 168 762 32 D-12*20.D0,
     9           b37=-0.420 024 676 982 08 D-05*20.D0,
     1           b38=-0.590 560 296 856 39 D-25*21.D0,
     2           b39= 0.378 269 476 134 57 D-05*22.D0,
     3           b40=-0.127 686 089 346 81 D-14*23.D0,
     4           b41= 0.730 876 105 950 61 D-28*24.D0,
     5           b42= 0.554 147 153 507 78 D-16*24.D0,
     6           b43=-0.943 697 072 412 10 D-06*24.D0)
      PARAMETER (a06=-0.330 326 416 702 03 D-04*2.D0,
     6           a07=-0.189 489 875 163 15 D-03*2.D0,
     7           a08=-0.393 927 772 433 55 D-02*2.D0,
     8           a09=-0.437 972 956 505 73 D-01*2.D0,
     9           a10=-0.266 745 479 140 87 D-04*2.D0,
     1           a11= 0.204 817 376 923 09 D-07*6.D0,
     2           a12= 0.438 706 672 844 35 D-06*6.D0,
     3           a13=-0.322 776 772 385 70 D-04*6.D0,
     4           a14=-0.150 339 245 421 48 D-02*6.D0,
     5           a15=-0.406 682 535 626 49 D-01*6.D0,
     6           a16=-0.788 473 095 593 67 D-09*12.D0,
     7           a17= 0.127 907 178 522 85 D-07*12.D0,
     8           a18= 0.482 253 727 185 07 D-06*12.D0,
     9           a19= 0.229 220 763 376 61 D-05*20.D0,
     1           a20=-0.167 147 664 510 61 D-10*30.D0,
     2           a21=-0.211 714 723 213 55 D-02*30.D0,
     3           a22=-0.238 957 419 341 04 D02*30.D0,
     4           a23=-0.590 595 643 242 70 D-17*42.D0,
     5           a24=-0.126 218 088 991 01 D-05*42.D0,
     6           a25=-0.389 468 424 357 39 D-01*42.D0,
     7           a26= 0.112 562 113 604 59 D-10*56.D0,
     8           a27=-0.823 113 408 979 98 D01*56.D0,
     9           a28= 0.198 097 128 020 88 D-07*72.D0,
     1           a29= 0.104 069 652 101 74 D-18*90.D0,
     2           a30=-0.102 347 470 959 29 D-12*90.D0,
     3           a31=-0.100 181 793 795 11 D-08*90.D0,
     4           a32=-0.808 829 086 469 85 D-10*240.D0,
     5           a33= 0.106 930 318 794 09 D0*240.D0,
     6           a34=-0.336 622 505 741 71 D0*306.D0,
     7           a35= 0.891 858 453 554 21 D-24*380.D0,
     8           a36= 0.306 293 168 762 32 D-12*380.D0,
     9           a37=-0.420 024 676 982 08 D-05*380.D0,
     1           a38=-0.590 560 296 856 39 D-25*420.D0,
     2           a39= 0.378 269 476 134 57 D-05*462.D0,
     3           a40=-0.127 686 089 346 81 D-14*506.D0,
     4           a41= 0.730 876 105 950 61 D-28*552.D0,
     5           a42= 0.554 147 153 507 78 D-16*552.D0,
     6           a43=-0.943 697 072 412 10 D-06*552.D0)
      PARAMETER (c03=-0.560 879 112 830 20 D-02*30.D0,
     3           c04= 0.714 527 380 814 55 D-01*20.D0,
     4           c05=-0.407 104 982 239 28 D0*12.D0,
     5           c06= 0.142 408 191 714 44 D01*6.D0,
     6           c07=-0.438 395 113 194 50 D01*2.D0,
     7           c08=-0.284 086 324 607 72 D0*2.D0,
     8           c09= 0.212 684 637 533 07 D-01*6.D0)
      PARAMETER (d03=-0.459 960 136 963 65 D-01*2.D0,
     3           d04=-0.575 812 590 834 32 D-01*6.D0,
     4           d05=-0.503 252 787 279 30 D-01*30.D0,
     6           d07=-0.189 489 875 163 15 D-03*2.D0,
     7           d08=-0.393 927 772 433 55 D-02*12.D0,
     8           d09=-0.437 972 956 505 73 D-01*42.D0,
     9           d10=-0.266 745 479 140 87 D-04*1260.D0,
     3           d13=-0.322 776 772 385 70 D-04*6.D0,
     4           d14=-0.150 339 245 421 48 D-02*30.D0,
     5           d15=-0.406 682 535 626 49 D-01*1190.D0,
     7           d17= 0.127 907 178 522 85 D-07*2.D0,
     8           d18= 0.482 253 727 185 07 D-06*6.D0,
     9           d19= 0.229 220 763 376 61 D-05*42.D0,
     1           d20=-0.167 147 664 510 61 D-10*6.D0,
     2           d21=-0.211 714 723 213 55 D-02*240.D0,
     3           d22=-0.238 957 419 341 04 D02*1190.D0,
     5           d24=-0.126 218 088 991 01 D-05*110.D0,
     6           d25=-0.389 468 424 357 39 D-01*600.D0,
     7           d26= 0.112 562 113 604 59 D-10*56.D0,
     8           d27=-0.823 113 408 979 98 D01*1260.D0,
     9           d28= 0.198 097 128 020 88 D-07*156.D0,
     1           d29= 0.104 069 652 101 74 D-18*12.D0,
     2           d30=-0.102 347 470 959 29 D-12*90.D0,
     3           d31=-0.100 181 793 795 11 D-08*182.D0,
     4           d32=-0.808 829 086 469 85 D-10*812.D0,
     5           d33= 0.106 930 318 794 09 D0*2450.D0,
     6           d34=-0.336 622 505 741 71 D0*3192.D0,
     7           d35= 0.891 858 453 554 21 D-24*380.D0,
     8           d36= 0.306 293 168 762 32 D-12*1190.D0,
     9           d37=-0.420 024 676 982 08 D-05*2256.D0,
     1           d38=-0.590 560 296 856 39 D-25*420.D0,
     2           d39= 0.378 269 476 134 57 D-05*2756.D0,
     3           d40=-0.127 686 089 346 81 D-14*1482.D0,
     4           d41= 0.730 876 105 950 61 D-28*650.D0,
     5           d42= 0.554 147 153 507 78 D-16*1560.D0,
     6           d43=-0.943 697 072 412 10 D-06*3306.D0)
      PARAMETER (e02=-0.178 348 622 923 58 D-01,
     2           e03=-0.459 960 136 963 65 D-01*2.D0,
     3           e04=-0.575 812 590 834 32 D-01*3.D0,
     4           e05=-0.503 252 787 279 30 D-01*6.D0,
     5           e06=-0.330 326 416 702 03 D-04*2.D0,
     6           e07=-0.189 489 875 163 15 D-03*4.D0,
     7           e08=-0.393 927 772 433 55 D-02*8.D0,
     8           e09=-0.437 972 956 505 73 D-01*14.D0,
     9           e10=-0.266 745 479 140 87 D-04*72.D0,
     2           e12= 0.438 706 672 844 35 D-06*3.D0,
     3           e13=-0.322 776 772 385 70 D-04*9.D0,
     4           e14=-0.150 339 245 421 48 D-02*18.D0,
     5           e15=-0.406 682 535 626 49 D-01*105.D0,
     6           e16=-0.788 473 095 593 67 D-09*4.D0,
     7           e17= 0.127 907 178 522 85 D-07*8.D0,
     8           e18= 0.482 253 727 185 07 D-06*12.D0,
     9           e19= 0.229 220 763 376 61 D-05*35.D0,
     1           e20=-0.167 147 664 510 61 D-10*18.D0,
     2           e21=-0.211 714 723 213 55 D-02*96.D0,
     3           e22=-0.238 957 419 341 04 D02*210.D0,
     5           e24=-0.126 218 088 991 01 D-05*77.D0,
     6           e25=-0.389 468 424 357 39 D-01*175.D0,
     7           e26= 0.112 562 113 604 59 D-10*64.D0,
     8           e27=-0.823 113 408 979 98 D01*288.D0,
     9           e28= 0.198 097 128 020 88 D-07*117.D0,
     1           e29= 0.104 069 652 101 74 D-18*40.D0,
     2           e30=-0.102 347 470 959 29 D-12*100.D0,
     3           e31=-0.100 181 793 795 11 D-08*140.D0,
     4           e32=-0.808 829 086 469 85 D-10*464.D0,
     5           e33= 0.106 930 318 794 09 D0*800.D0,
     6           e34=-0.336 622 505 741 71 D0*1026.D0,
     7           e35= 0.891 858 453 554 21 D-24*400.D0,
     8           e36= 0.306 293 168 762 32 D-12*700.D0,
     9           e37=-0.420 024 676 982 08 D-05*960.D0,
     1           e38=-0.590 560 296 856 39 D-25*441.D0,
     2           e39= 0.378 269 476 134 57 D-05*1166.D0,
     3           e40=-0.127 686 089 346 81 D-14*897.D0,
     4           e41= 0.730 876 105 950 61 D-28*624.D0,
     5           e42= 0.554 147 153 507 78 D-16*960.D0,
     6           e43=-0.943 697 072 412 10 D-06*1392.D0)
C
      pi= p
      tau=540.d0/t
      tauinv= 1.d0/tau
C
      taun= tau - 0.5d0
      taun2= taun*taun
      taun3= taun2*taun
      taun4= taun3*taun
      taun6= taun4*taun2
      taun7= taun6*taun
      taun9= taun7*taun2
      taun10= taun9*taun
      taun13= taun10*taun3
      taun14= taun13*taun
      taun15= taun14*taun
      taun19= taun15*taun4
      taun20= taun19*taun
      taun29= taun20*taun9
C
      pi2= pi*pi
      pi4= pi2*pi2
C
      gap2= b01 + taun*(b02 + taun*(b03 + taun*(b04 + taun3*b05))) +
     3     pi* (taun*(b06 + taun*(b07 + taun2*(b08 +
     4         taun3*(b09 + taun29*b10)))) +
     5     pi* (b11 + taun*(b12 + taun2*(b13 + taun3*(b14 +
     6         taun29*b15))) +
     7     pi* (taun*(b16 + taun*(b17 + taun*b18)) +
     8     pi* (taun7*b19 +
     9     pi* (taun3*(b20 + taun13*(b21 + taun19*b22)) +
     1     pi* (b23 + taun10*taun*(b24 + taun14*b25) +
     2     pi* (taun7*(taun*b26 + taun29*b27) +
     3     pi* (taun13*b28 +
     4     pi* (taun4*(b29 + taun6*(b30 + taun4*b31)) +
     5     pi4*pi2*taun20*(taun9*b32 + taun20*taun10*(b33 +
     6         pi2*taun7*b34) +
     7     pi4*(b35 + taun15*(b36 + taun13*b37) +
     8     pi* (taun*b38 +
     9     pi* (taun19*(taun14*b39 +
     1          pi* b40) +
     2     pi2*(taun6*(b41 + taun14*(b42 +
     3     taun14*taun4*b43))) ))) ))) ))) ))) )
C
      gapp2= (taun*(a06 + taun*(a07 + taun2*(a08 +
     4     taun3*(a09 + taun29*a10)))) +
     5     pi* (a11 + taun*(a12 + taun2*(a13 + taun3*(a14 +
     6         taun29*a15))) +
     7     pi* (taun*(a16 + taun*(a17 + taun*a18)) +
     8     pi* (taun7*a19 +
     9     pi* (taun3*(a20 + taun13*(a21 + taun19*a22)) +
     1     pi* (a23 + taun10*taun*(a24 + taun14*a25) +
     2     pi* (taun7*(taun*a26 + taun29*a27) +
     3     pi* (taun13*a28 +
     4     pi* (taun4*(a29 + taun6*(a30 + taun4*a31)) +
     5     pi4*pi2*taun20*(taun9*a32 + taun20*taun10*(a33 +
     6         pi2*taun7*a34) +
     7     pi4*(a35 + taun15*(a36 + taun13*a37) +
     8     pi* (taun*a38 +
     9     pi* (taun19*(taun14*a39 +
     1          pi* a40) +
     2     pi2*(taun6*(a41 + taun14*(a42 +
     3     taun14*taun4*a43))) ))) ))) ))) ))) )
C
      gatt2= c08 + tau*c09 +
     1     tauinv*tauinv*tauinv*(c07 + tauinv*(c06 +
     2     tauinv*(c05 + tauinv*(c04 + tauinv*c03)))) +
     3     pi* (d03 + taun*(d04 + taun3*d05) +
     3     pi* (d07 + taun2*(d08 + taun3*(d09 + taun29*d10)) +
     5     pi* (taun*(d13 + taun3*(d14 + taun29*d15)) +
     7     pi* (d17 + taun*d18 +
     8     pi* (taun3*taun2*d19 +
     9     pi* (taun*(d20 + taun13*(d21 + taun19*d22)) +
     1     pi* (taun9*(d24 + taun14*d25) +
     2     pi* (taun3*taun2*(taun*d26 + taun29*d27) +
     3     pi* (taun10*taun*d28 +
     4     pi* (taun2*(d29 + taun6*(d30 + taun4*d31)) +
     5     pi4*pi2*taun15*taun3*(taun9*d32 + taun20*taun10*(d33 +
     6         pi2*taun7*d34) +
     7     pi4*(d35 + taun15*(d36 + taun13*d37) +
     8     pi* (taun*d38 +
     9     pi* (taun19*(taun14*d39 + pi* d40) +
     2     pi2*(taun6*(d41 + taun14*(d42 +
     3     taun14*taun4*d43))) ))) ))) ))) ))) ))
C
      gapt2= e02 + taun*(e03 + taun*(e04 + taun3*e05)) +
     3     pi* (e06 + taun*(e07 + taun2*(e08 +
     4         taun3*(e09 + taun29*e10))) +
     5     pi* (e12 + taun2*(e13 + taun3*(e14 + taun29*e15)) +
     7     pi* (e16 + taun*(e17 + taun*e18) +
     8     pi* (taun6*e19 +
     9     pi* (taun2*(e20 + taun13*(e21 + taun19*e22)) +
     1     pi* (taun10*(e24 + taun14*e25) +
     2     pi* (taun6*(taun*e26 + taun29*e27) +
     3     pi* (taun6*taun6*e28 +
     4     pi* (taun3*(e29 + taun6*(e30 + taun4*e31)) +
     5     pi4*pi2*taun19*(taun9*e32 + taun20*taun10*(e33 +
     6         pi2*taun7*e34) +
     7     pi4*(e35 + taun15*(e36 + taun13*e37) +
     8     pi* (taun*e38 +
     9     pi* (taun19*(taun14*e39 + pi* e40) +
     2     pi2*(taun6*(e41 + taun14*(e42 +
     3     taun14*taun4*e43))) ))) ))) ))) ))) )
C
      cvpt2n= 0.461526d0*(-tau*tau*gatt2-(1.d0+pi*gap2-tau*pi*gapt2)*
     1           (1.d0+pi*gap2-tau*pi*gapt2)/(1.d0-pi2*gapp2))
C
      END
C
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is specific volume in m^3/kg.
C  Second argument is temperature in K.
C  Returns specific enthalpy in kJ/kg
c
      FUNCTION HVT3N(V,T)
c
      IMPLICIT double precision (A-H,O-Z)
c
      parameter (f1= 0.49189764279794D+00,
     &   f3 =  0.96663837579639D+01,
     &   f4 = -0.70952891492865D+01,
     &   f5 =  0.84598470171481D+01,
     &   f6 = -0.12960010600398D+02,
     &   f7 =  0.66755322030656D+01,
     &   f8 = -0.89768500555238D-01,
     &   f9 = -0.17520886815502D+01,
     &  f10 = -0.37231696861705D+01,
     &  f11 =  0.65367621353450D+01,
     &  f12 = -0.53340395459772D+01,
     &  f13 =  0.35531465412227D+00,
     &  f14 = -0.15731521481920D+01,
     &  f15 =  0.18081584968717D+02,
     &  f16 = -0.12669975838924D+02,
     &  f17 =  0.43664646488852D+00,
     &  f18 =  0.16228889618024D+01,
     &  f19 = -0.38767255915581D+00,
     &  f20 =  0.32075594480472D+01,
     &  f21 = -0.65227195048178D+01,
     &  f22 = -0.72035213749267D-01,
     &  f23 = -0.63703743166930D+01,
     &  f24 =  0.81199175821825D-01,
     &  f25 = -0.12316218802611D+01,
     &  f26 =  0.33441092408796D+01,
     &  f27 =  0.97643832897896D+01,
     &  f28 =  0.29825229662283D+00,
     &  f29 = -0.12152394310570D+01,
     &  f30 = -0.72782785116784D+01,
     &  f31 = -0.61407144380140D-01,
     &  f32 =  0.34803030138386D+00,
     &  f33 =  0.24274463510316D+01,
     &  f34 = -0.56089416393544D-01,
     &  f35 = -0.23277882648693D+00,
     &  f36 =  0.29406244005770D-02,
     &  f37 =  0.52189898432996D-01,
     &  f38 =  0.37367361667631D-03,
     &  f39 = -0.84059796975928D-03,
     &  f40 = -0.76714125522093D-03)
C
      dn=1.d0/(v*322.d0)
      dn2=dn*dn
c
      tn=647.096d0/t
      tn2=tn*tn
      tn4=tn2*tn2
      tn5=tn4*tn
      tn8=tn4*tn4
      tn12=tn8*tn4
C
      HVT3N = ((((((((((((((dn2*f40 + f37)*dn + f35)*dn2 + f33)*dn +
     &   f30)*dn + f27)*dn + f23)*dn + f18)*tn4 + f17)*dn*tn5 +
     &   f12)*tn2 + f11)*tn*tn8 + f10)*tn4 + f9)*tn2 + ((((((((tn*f39
     &   + f38)*dn + tn2*f36)*dn2 + tn2*f34)*dn + tn2*f32 + f31)*
     &   dn + (tn2*f29 + f28)*tn)*dn + (tn2*f26 + f25)*
     &   tn2 + f24)*dn + ((tn12*f22 + f21)*tn2 + f20)*tn2 + f19)*
     &   dn + ((tn*f16 + f15)*tn4 + f14)*tn2 + f13)*dn)*dn +
     &   (((tn*tn12*f8 + tn2*f7 + f6)*tn8 + tn5*f5 + f4)*tn + f3)*tn +
     &   f1)*t
c
      END
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is specific volume in m^3/kg.
C  Second argument is temperature in K.
C  Returns pressure in MPa.
c
      FUNCTION PVT3N(V,T)
c
      IMPLICIT double precision (A-H,O-Z)
c
      parameter (f1= 0.49189764279794D-03,
     &   f9 = -0.18137563991203D-05,
     &  f10 = -0.16518055395610D-05,
     &  f11 =  0.12687814703698D-05,
     &  f12 = -0.92029667804989D-06,
     &  f13 =  0.11034616587648D-05,
     &  f14 = -0.24427828388075D-05,
     &  f15 =  0.14038497646519D-04,
     &  f16 = -0.87439446783466D-05,
     &  f17 =  0.11300374350117D-06,
     &  f18 =  0.36000198797746D-06,
     &  f19 = -0.12039520470677D-05,
     &  f20 =  0.59768188472930D-05,
     &  f21 = -0.86815255166165D-05,
     &  f22 = -0.35322922727656D-07,
     &  f23 = -0.20465970175711D-05,
     &  f24 =  0.25217135348393D-06,
     &  f25 = -0.25499417810788D-05,
     &  f26 =  0.51927162125460D-05,
     &  f27 =  0.40432228943229D-05,
     &  f28 =  0.77187447366156D-06,
     &  f29 = -0.23587721876106D-05,
     &  f30 = -0.36457015185726D-05,
     &  f31 = -0.19070541732963D-06,
     &  f32 =  0.81062958396863D-06,
     &  f33 =  0.14134974870137D-05,
     &  f34 = -0.13548168211001D-06,
     &  f35 = -0.17009779063714D-06,
     &  f36 =  0.74719422939563D-08,
     &  f37 =  0.41677824835579D-07,
     &  f38 =  0.11604770704233D-08,
     &  f39 = -0.23732297282870D-08,
     &  f40 = -0.70828888764732D-09)
C
      d=1.d0/v
      dn=d*0.31055900621118d-2
      dn2=dn*dn
c
      tn=647.096d0/t
      tn2=tn*tn
      tn4=tn2*tn2
      tn8=tn4*tn4
C
      PVT3N = (((((((((((((dn2*f40 + f37)*dn + f35)*dn2 + f33)*dn +
     &   f30)*dn + f27)*dn + f23)*dn + f18)*tn4 + f17)*dn*tn*tn4 +
     &   f12)*tn2 + f11)*tn*tn8 + f10)*tn4 + f9)*tn2 + ((((((((tn*
     &   f39 + f38)*dn + tn2*f36)*dn2 + tn2*f34)*dn + tn2*f32 + f31)*
     &   dn + (tn2*f29 + f28)*tn)*dn + (tn2*f26 + f25)*
     &   tn2 + f24)*dn + ((tn4*tn8*f22 + f21)*tn2 + f20)*tn2 + f19)*
     &   dn + ((tn*f16 + f15)*tn4 + f14)*tn2 + f13)*dn + v*f1)*t*d*d
c
      END
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is specific volume in m^3/kg.
C  Second argument is temperature in K.
C  Returns specific entropy in kJ/(kg*K)
c
      FUNCTION SVT3N(V,T)
c
      IMPLICIT double precision (A-H,O-Z)
c
      parameter (f1= -0.49189764279795D+00,
     &   f2 =  0.72611171554229D+01,
     &   f4 = -0.35476445746432D+01,
     &   f5 =  0.72512974432699D+01,
     &   f6 = -0.11664009540359D+02,
     &   f7 =  0.61192378528101D+01,
     &   f8 = -0.85865522270227D-01,
     &   f9 = -0.58402956051675D+00,
     &  f10 = -0.26594069186932D+01,
     &  f11 =  0.57196668684268D+01,
     &  f12 = -0.47413684853130D+01,
     &  f13 = -0.17765732706114D+00,
     &  f14 = -0.39328803704801D+00,
     &  f15 =  0.11300990605448D+02,
     &  f16 = -0.84466505592829D+01,
     &  f17 =  0.38206565677745D+00,
     &  f18 =  0.14490080016093D+01,
     &  f19 =  0.12922418638527D+00,
     &  f20 =  0.64151188960945D+00,
     &  f21 = -0.27954512163505D+01,
     &  f22 = -0.56869905591526D-01,
     &  f23 = -0.54917019971491D+01,
     &  f24 = -0.20299793955457D-01,
     &  f25 = -0.20527031337684D+00,
     &  f26 =  0.12540409653298D+01,
     &  f27 =  0.81369860748247D+01,
     &  f29 = -0.30380985776425D+00,
     &  f30 = -0.58695794449019D+01,
     &  f31 =  0.10234524063357D-01,
     &  f32 =  0.43503787672983D-01,
     &  f33 =  0.18964424617434D+01,
     &  f34 = -0.62321573770604D-02,
     &  f35 = -0.17116090182863D+00,
     &  f36 =  0.26732949096155D-03,
     &  f37 =  0.37278498880712D-01,
     &  f38 = -0.37367361667631D-04,
     &  f40 = -0.51833868596009D-03)
C
      dn=1.d0/(v*322.d0)
      dn2=dn*dn
c
      tn=647.096d0/t
      tn2=tn*tn
      tn4=tn2*tn2
      tn5=tn4*tn
      tn8=tn4*tn4
      tn12=tn8*tn4
C
      SVT3N = (((((((((((((dn2*f40 + f37)*dn + f35)*dn2 + f33)*dn +
     &   f30)*dn + f27)*dn + f23)*dn + f18)*tn4 + f17)*dn*tn5 +
     &   f12)*tn2 + f11)*tn*tn8 + f10)*tn4 + f9)*tn2 + (((((((
     &   f38*dn + tn2*f36)*dn2 + tn2*f34)*dn + tn2*f32 + f31)*
     &   dn + tn2*f29*tn)*dn + (tn2*f26 + f25)*
     &   tn2 + f24)*dn + ((tn12*f22 + f21)*tn2 + f20)*tn2 + f19)*
     &   dn + ((tn*f16 + f15)*tn4 + f14)*tn2 + f13)*dn)*dn +
     &   ((tn*tn12*f8 + tn2*f7 + f6)*tn8 + tn5*f5 + f4)*tn2 +
     &   f1*dlog(dn) + f2
c
      END
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is specific volume in m^3/kg.
C  Second argument is temperature in K.
C  Returns speed of sound in m/s
C
      FUNCTION WVT3N(V,T)
C
      IMPLICIT double precision (A-H,O-Z)
C
      PARAMETER (a01=  0.106 580 700 285 13 D01,
     1           a09= -0.126 543 154 777 14 D01,
     9           a10= -0.115 244 078 066 81 D01,
     1           a11=  0.885 210 439 843 18 D0,
     2           a12= -0.642 077 651 816 07 D0,
     3           a13=  0.384 934 601 866 71 D0*2.D0,
     4           a14= -0.852 147 088 242 06 D0*2.D0,
     5           a15=  0.489 722 815 418 77 D01*2.D0,
     6           a16= -0.305 026 172 569 65 D01*2.D0,
     7           a17=  0.394 205 368 791 54 D-01*2.D0,
     8           a18=  0.125 584 084 243 08 D0*2.D0,
     9           a19= -0.279 993 296 987 10 D0*3.D0,
     1           a20=  0.138 997 995 694 60 D01*3.D0,
     2           a21= -0.201 899 150 235 70 D01*3.D0,
     3           a22= -0.821 476 371 739 63 D-02*3.D0,
     4           a23= -0.475 960 357 349 23 D0*3.D0,
     5           a24=  0.439 840 744 735 00 D-01*4.D0,
     6           a25= -0.444 764 354 287 39 D0*4.D0,
     7           a26=  0.905 720 707 197 33 D0*4.D0,
     8           a27=  0.705 224 500 879 67 D0*4.D0,
     9           a28=  0.107 705 126 263 32 D0*5.D0,
     1           a29= -0.329 136 232 589 54 D0*5.D0,
     2           a30= -0.508 710 620 411 58 D0*5.D0,
     3           a31= -0.221 754 008 730 96 D-01*6.D0,
     4           a32=  0.942 607 516 650 92 D-01*6.D0,
     5           a33=  0.164 362 784 479 61 D0*6.D0,
     6           a34= -0.135 033 722 413 48 D-01*7.D0,
     7           a35= -0.148 343 453 524 72 D-01*8.D0,
     8           a36=  0.579 229 536 280 84 D-03*9.D0,
     9           a37=  0.323 089 047 037 11 D-02*9.D0,
     1           a38=  0.809 648 029 962 15 D-04*10.D0,
     2           a39= -0.165 576 797 950 37 D-03*10.D0,
     3           a40= -0.449 238 990 618 15 D-04*11.D0)
      PARAMETER (b01=  0.106 580 700 285 13 D01,
     3           b13=  0.384 934 601 866 71 D0*2.D0,
     4           b14= -0.852 147 088 242 06 D0*2.D0,
     5           b15=  0.489 722 815 418 77 D01*2.D0,
     6           b16= -0.305 026 172 569 65 D01*2.D0,
     7           b17=  0.394 205 368 791 54 D-01*2.D0,
     8           b18=  0.125 584 084 243 08 D0*2.D0,
     9           b19= -0.279 993 296 987 10 D0*6.D0,
     1           b20=  0.138 997 995 694 60 D01*6.D0,
     2           b21= -0.201 899 150 235 70 D01*6.D0,
     3           b22= -0.821 476 371 739 63 D-02*6.D0,
     4           b23= -0.475 960 357 349 23 D0*6.D0,
     5           b24=  0.439 840 744 735 00 D-01*12.D0,
     6           b25= -0.444 764 354 287 39 D0*12.D0,
     7           b26=  0.905 720 707 197 33 D0*12.D0,
     8           b27=  0.705 224 500 879 67 D0*12.D0,
     9           b28=  0.107 705 126 263 32 D0*20.D0,
     1           b29= -0.329 136 232 589 54 D0*20.D0,
     2           b30= -0.508 710 620 411 58 D0*20.D0,
     3           b31= -0.221 754 008 730 96 D-01*30.D0,
     4           b32=  0.942 607 516 650 92 D-01*30.D0,
     5           b33=  0.164 362 784 479 61 D0*30.D0,
     6           b34= -0.135 033 722 413 48 D-01*42.D0,
     7           b35= -0.148 343 453 524 72 D-01*56.D0,
     8           b36=  0.579 229 536 280 84 D-03*72.D0,
     9           b37=  0.323 089 047 037 11 D-02*72.D0,
     1           b38=  0.809 648 029 962 15 D-04*90.D0,
     2           b39= -0.165 576 797 950 37 D-03*90.D0,
     3           b40= -0.449 238 990 618 15 D-04*110.D0)
      PARAMETER (c04= -0.768 677 078 787 16 D01*2.D0,
     4           c05=  0.261 859 477 879 54 D01*42.D0,
     5           c06= -0.280 807 811 486 20 D01*90.D0,
     6           c07=  0.120 533 696 965 17 D01*132.D0,
     7           c08= -0.845 668 128 125 02 D-02*506.D0,
     8           c09= -0.126 543 154 777 14 D01*2.D0,
     9           c10= -0.115 244 078 066 81 D01*30.D0,
     1           c11=  0.885 210 439 843 18 D0*210.D0,
     2           c12= -0.642 077 651 816 07 D0*272.D0,
     4           c14= -0.852 147 088 242 06 D0*2.D0,
     5           c15=  0.489 722 815 418 77 D01*30.D0,
     6           c16= -0.305 026 172 569 65 D01*42.D0,
     7           c17=  0.394 205 368 791 54 D-01*462.D0,
     8           c18=  0.125 584 084 243 08 D0*650.D0,
     1           c20=  0.138 997 995 694 60 D01*2.D0,
     2           c21= -0.201 899 150 235 70 D01*12.D0,
     3           c22= -0.821 476 371 739 63 D-02*240.D0,
     4           c23= -0.475 960 357 349 23 D0*650.D0,
     6           c25= -0.444 764 354 287 39 D0*2.D0,
     7           c26=  0.905 720 707 197 33 D0*12.D0,
     8           c27=  0.705 224 500 879 67 D0*650.D0,
     1           c29= -0.329 136 232 589 54 D0*6.D0,
     2           c30= -0.508 710 620 411 58 D0*650.D0,
     4           c32=  0.942 607 516 650 92 D-01*2.D0,
     5           c33=  0.164 362 784 479 61 D0*650.D0,
     6           c34= -0.135 033 722 413 48 D-01*2.D0,
     7           c35= -0.148 343 453 524 72 D-01*650.D0,
     8           c36=  0.579 229 536 280 84 D-03*2.D0,
     9           c37=  0.323 089 047 037 11 D-02*650.D0,
     3           c40= -0.449 238 990 618 15 D-04*650.D0)
      PARAMETER (d09= -0.126 543 154 777 14 D01*2.D0,
     9           d10= -0.115 244 078 066 81 D01*6.D0,
     1           d11=  0.885 210 439 843 18 D0*15.D0,
     2           d12= -0.642 077 651 816 07 D0*17.D0,
     4           d14= -0.852 147 088 242 06 D0*4.D0,
     5           d15=  0.489 722 815 418 77 D01*12.D0,
     6           d16= -0.305 026 172 569 65 D01*14.D0,
     7           d17=  0.394 205 368 791 54 D-01*44.D0,
     8           d18=  0.125 584 084 243 08 D0*52.D0,
     1           d20=  0.138 997 995 694 60 D01*6.D0,
     2           d21= -0.201 899 150 235 70 D01*12.D0,
     3           d22= -0.821 476 371 739 63 D-02*48.D0,
     4           d23= -0.475 960 357 349 23 D0*78.D0,
     6           d25= -0.444 764 354 287 39 D0*8.D0,
     7           d26=  0.905 720 707 197 33 D0*16.D0,
     8           d27=  0.705 224 500 879 67 D0*104.D0,
     9           d28=  0.107 705 126 263 32 D0*5.D0,
     1           d29= -0.329 136 232 589 54 D0*15.D0,
     2           d30= -0.508 710 620 411 58 D0*130.D0,
     4           d32=  0.942 607 516 650 92 D-01*12.D0,
     5           d33=  0.164 362 784 479 61 D0*156.D0,
     6           d34= -0.135 033 722 413 48 D-01*14.D0,
     7           d35= -0.148 343 453 524 72 D-01*208.D0,
     8           d36=  0.579 229 536 280 84 D-03*18.D0,
     9           d37=  0.323 089 047 037 11 D-02*234.D0,
     2           d39= -0.165 576 797 950 37 D-03*10.D0,
     3           d40= -0.449 238 990 618 15 D-04*286.D0)
C
      PARAMETER ( tn= 647.096D0,
     2            rhonq= 1.D0/322.0D0,
     3            r= 0.461526D0,
     4            rnn= r*1000.0D0) 
C
      rho= 1.d0/v
      tau= tn/t 
      del= rho*rhonq
C
      tau2= tau*tau
      tau4= tau2*tau2
      tau12= tau4*tau4*tau4
C
      del2= del*del
C
      phd3= a01/del + tau2*(a09 + tau4*a10 + tau12*tau*(a11 +
     1      tau2*(a12 + del*tau4*tau*(a17 + tau4*(a18 + del*(a23 +
     2      del*(a27 + del*(a30 + del*(a33 + del2*(a35 +
     3      del*(a37 + del2*a40))) ))) ))) )) +
     4      del* (a13 + tau2*(a14 + tau4*(a15 + tau*a16)) +
     6      del* (a19 + tau2*(a20 + tau2*(a21 + tau12*a22)) +
     7      del* (a24 + tau2*(a25 + tau2*a26) +
     8      del* (tau*(a28 + tau2*a29) +
     1      del* (a31 + tau2*a32 +
     2      del* (tau2*a34 +
     3      del2*(tau2*a36 +
     4      del* (a38 + tau*a39))) ))) ))
C
      phdd3= - b01/del2 + b13 + tau2*(b14 + tau4*(b15 + tau*b16 +
     1       tau12*tau4*(b17 + tau4*(b18 + del*(b23 + del*(b27 +
     2       del*(b30 + del*(b33 + del2*(b35 + del*(b37 +
     3       del2*b40))) ))) ))) ) +
     6       del* (b19 + tau2*(b20 + tau2*(b21 + tau12*b22)) +
     7       del* (b24 + tau2*(b25 + tau2*b26) +
     8       del* (tau*(b28 + tau2*b29) +
     1       del* (b31 + tau2*b32 +
     2       del* (tau2*b34 +
     3       del2*(tau2*b36 +
     4       del* (b38 + tau*b39))) ))) )
C
      phtt3= c04 + tau4*(tau*c05 + tau4*(c06 + tau2*c07 +
     1       tau12*(tau*c08 + del2*(c17 + tau4*(c18 + del*(c23 +
     2       del*(c27 + del*(c30 + del*(c33 + del2*(c35 +
     4       del*(c37 + del2*c40))) ))) ))) )) +
     2       del* (c09 + tau4*c10 + tau12*tau*(c11 + tau2*c12) +
     4       del* (c14 + tau4*(c15 + tau*c16) +
     6       del* (c20 + tau2*(c21 + tau12*c22) +
     7       del* (c25 + tau2*c26 +
     8       del* (tau*c29 +
     1       del* (c32 +
     2       del* (c34 +
     3       del2*c36))) ))) )
C
      phdt3= tau*(d09 + tau4*d10 + tau12*tau*(d11 +
     1       tau2*(d12 + del*tau4*tau*(d17 + tau4*(d18 + del*(d23 +
     2       del*(d27 + del*(d30 + del*(d33 + del2*(d35 +
     3       del*(d37 + del2*d40))) ))) ))) )) +
     4       del* (tau*(d14 + tau4*(d15 + tau*d16)) +
     6       del* (tau*(d20 + tau2*(d21 + tau12*d22)) +
     7       del* (tau*(d25 + tau2*d26) +
     8       del* (d28 + tau2*d29 +
     1       del* (tau*d32 +
     2       del* (tau*d34 +
     3       del2*(tau*d36 +
     4       del*d39))) ))) )
C      
      wvt3n= DSQRT(rnn*t*(2.d0*del*phd3 + del2*phdd3 -
     1       (del*phd3 - del*tau*phdt3)*(del*phd3 - del*tau*phdt3)/
     2       (tau2*phtt3)))
C
      END
C
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is specific volume in m^3/kg.
C  Second argument is temperature in K.
C  Returns specific isobaric heat capacity in kJ/(kg*K)
c
      FUNCTION CPVT3N(V,T)
c
      IMPLICIT double precision (A-H,O-Z)
c
      parameter (e1= 0.10658070028513D+04,
     &   e9 =  0.12654315477714D+04,
     &  e10 =  0.57622039033404D+04,
     &  e11 = -0.12392946157805D+05,
     &  e12 =  0.10273242429057D+05,
     &  e13 =  0.76986920373342D+03,
     &  e14 =  0.17042941764841D+04,
     &  e15 = -0.48972281541877D+05,
     &  e16 =  0.36603140708358D+05,
     &  e17 = -0.16556625489245D+04,
     &  e18 = -0.62792042121540D+04,
     &  e19 = -0.83997989096131D+03,
     &  e20 = -0.41699398708380D+04,
     &  e21 =  0.18170923521213D+05,
     &  e22 =  0.36966436728283D+03,
     &  e23 =  0.35697026801192D+05,
     &  e24 =  0.17593629789400D+03,
     &  e25 =  0.17790574171496D+04,
     &  e26 = -0.10868648486368D+05,
     &  e27 = -0.70522450087968D+05,
     &  e29 =  0.32913623258954D+04,
     &  e30 =  0.63588827551446D+05,
     &  e31 = -0.13305240523858D+03,
     &  e32 = -0.56556450999056D+03,
     &  e33 = -0.24654417671942D+05,
     &  e34 =  0.94523605689435D+02,
     &  e35 =  0.29668690704943D+04,
     &  e36 = -0.52130658265276D+01,
     &  e37 = -0.72695035583351D+03,
     &  e38 =  0.80964802996214D+00,
     &  e40 =  0.12354072241999D+02)
c
      parameter (f1= 0.10658070028513D+04,
     &   f9 = -0.25308630955428D+04,
     &  f10 = -0.23048815613362D+04,
     &  f11 =  0.17704208796864D+04,
     &  f12 = -0.12841553036321D+04,
     &  f13 =  0.23096076112002D+04,
     &  f14 = -0.51128825294524D+04,
     &  f15 =  0.29383368925126D+05,
     &  f16 = -0.18301570354179D+05,
     &  f17 =  0.23652322127492D+03,
     &  f18 =  0.75350450545848D+03,
     &  f19 = -0.33599195638452D+04,
     &  f20 =  0.16679759483352D+05,
     &  f21 = -0.24227898028284D+05,
     &  f22 = -0.98577164608756D+02,
     &  f23 = -0.57115242881908D+04,
     &  f24 =  0.87968148946999D+03,
     &  f25 = -0.88952870857479D+04,
     &  f26 =  0.18114414143946D+05,
     &  f27 =  0.14104490017593D+05,
     &  f28 =  0.32311537878995D+04,
     &  f29 = -0.98740869776863D+04,
     &  f30 = -0.15261318612347D+05,
     &  f31 = -0.93136683667005D+03,
     &  f32 =  0.39589515699338D+04,
     &  f33 =  0.69032369481436D+04,
     &  f34 = -0.75618884551549D+03,
     &  f35 = -0.10680728653780D+04,
     &  f36 =  0.52130658265275D+02,
     &  f37 =  0.29078014233339D+03,
     &  f38 =  0.89061283295836D+01,
     &  f39 = -0.18213447774541D+02,
     &  f40 = -0.59299546761597D+01)
C
      parameter ( g4 = 0.15373541575743D+05,
     &   g5 = -0.10998098070941D+06,
     &   g6 =  0.25272703033758D+06,
     &   g7 = -0.15910447999402D+06,
     &   g8 =  0.42790807283126D+04,
     &   g9 =  0.25308630955428D+04,
     &  g10 =  0.34573223420043D+05,
     &  g11 = -0.18589419236707D+06,
     &  g12 =  0.17464512129397D+06,
     &  g14 =  0.17042941764841D+04,
     &  g15 = -0.14691684462563D+06,
     &  g16 =  0.12811099247925D+06,
     &  g17 = -0.18212288038169D+05,
     &  g18 = -0.81629654758002D+05,
     &  g20 = -0.27799599138920D+04,
     &  g21 =  0.24227898028284D+05,
     &  g22 =  0.19715432921751D+04,
     &  g23 =  0.30937423227700D+06,
     &  g25 =  0.88952870857478D+03,
     &  g26 = -0.10868648486368D+05,
     &  g27 = -0.45839592557178D+06,
     &  g29 =  0.19748173955372D+04,
     &  g30 =  0.33066190326753D+06,
     &  g32 = -0.18852150333018D+03,
     &  g33 = -0.10683580991175D+06,
     &  g34 =  0.27006744482696D+02,
     &  g35 =  0.96423244791068D+04,
     &  g36 = -0.11584590725617D+01,
     &  g37 = -0.21000788057412D+04,
     &  g40 =  0.29200534390180D+02)
C
      parameter(R=0.461526d-3)
C
      dn=1.d0/(v*322.d0)
      dn2=dn*dn
c
      tn=647.096d0/t
      tn2=tn*tn
      tn3=tn2*tn
      tn4=tn3*tn
      tn5=tn4*tn
      tn8=tn4*tn4
      tn9=tn8*tn
      tn12=tn9*tn3
C
      FDN = (((((((((((((dn2*f40 + f37)*dn + f35)*dn2 + f33)*dn +
     &   f30)*dn + f27)*dn + f23)*dn + f18)*tn4 + f17)*dn*tn5 +
     &   f12)*tn2 + f11)*tn9 + f10)*tn4 + f9)*tn2 + ((((((((tn*
     &   f39 + f38)*dn + tn2*f36)*dn2 + tn2*f34)*dn + tn2*f32 + f31)*
     &   dn + (tn2*f29 + f28)*tn)*dn + (tn2*f26 + f25)*
     &   tn2 + f24)*dn + ((tn12*f22 + f21)*tn2 + f20)*tn2 + f19)*
     &   dn + ((tn*f16 + f15)*tn4 + f14)*tn2 + f13)*dn)*dn + f1
c
      FTTN = (((((((((((((dn2*g40 + g37)*dn + g35)*dn2 + g33)*dn +
     &   g30)*dn + g27)*dn + g23)*dn + g18)*tn4 + g17)*dn*tn5 +
     &   g12)*tn2 + g11)*tn9 + g10)*tn4 + g9)*tn2 + ((((((
     &   tn2*g36*dn2 + tn2*g34)*dn + tn2*g32)*
     &   dn + tn3*g29)*dn + (tn2*g26 + g25)*
     &   tn2)*dn + ((tn12*g22 + g21)*tn2 + g20)*tn2)*
     &   dn + ((tn*g16 + g15)*tn4 + g14)*tn2)*dn)*dn +
     &   ((tn*tn12*g8 + tn2*g7 + g6)*tn8 + tn5*g5 + g4)*tn2
c
      FDTN = (((((((((((((dn2*e40 + e37)*dn + e35)*dn2 + e33)*dn +
     &   e30)*dn + e27)*dn + e23)*dn + e18)*tn4 + e17)*dn*tn5 +
     &   e12)*tn2 + e11)*tn9 + e10)*tn4 + e9)*tn2 + (((((((
     &   e38*dn + tn2*e36)*dn2 + tn2*e34)*dn + tn2*e32 + e31)*
     &   dn + tn3*e29)*dn + (tn2*e26 + e25)*
     &   tn2 + e24)*dn + ((tn12*e22 + e21)*tn2 + e20)*tn2 + e19)*
     &   dn + ((tn*e16 + e15)*tn4 + e14)*tn2 + e13)*dn)*dn + e1
c
      CPVT3N = (FDTN/FDN*FDTN + FTTN) * R
c
      END
c
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is specific volume in m^3/kg.
C  Second argument is temperature in K.
C  Returns specific Gibbs free energy in kJ/kg
C
      FUNCTION GVT3N(V,T)
C
      IMPLICIT double precision (A-H,O-Z)
C
      PARAMETER (a01=  0.106 580 700 285 13 D01,
     1           a09= -0.126 543 154 777 14 D01,
     9           a10= -0.115 244 078 066 81 D01,
     1           a11=  0.885 210 439 843 18 D0,
     2           a12= -0.642 077 651 816 07 D0,
     3           a13=  0.384 934 601 866 71 D0*2.D0,
     4           a14= -0.852 147 088 242 06 D0*2.D0,
     5           a15=  0.489 722 815 418 77 D01*2.D0,
     6           a16= -0.305 026 172 569 65 D01*2.D0,
     7           a17=  0.394 205 368 791 54 D-01*2.D0,
     8           a18=  0.125 584 084 243 08 D0*2.D0,
     9           a19= -0.279 993 296 987 10 D0*3.D0,
     1           a20=  0.138 997 995 694 60 D01*3.D0,
     2           a21= -0.201 899 150 235 70 D01*3.D0,
     3           a22= -0.821 476 371 739 63 D-02*3.D0,
     4           a23= -0.475 960 357 349 23 D0*3.D0,
     5           a24=  0.439 840 744 735 00 D-01*4.D0,
     6           a25= -0.444 764 354 287 39 D0*4.D0,
     7           a26=  0.905 720 707 197 33 D0*4.D0,
     8           a27=  0.705 224 500 879 67 D0*4.D0,
     9           a28=  0.107 705 126 263 32 D0*5.D0,
     1           a29= -0.329 136 232 589 54 D0*5.D0,
     2           a30= -0.508 710 620 411 58 D0*5.D0,
     3           a31= -0.221 754 008 730 96 D-01*6.D0,
     4           a32=  0.942 607 516 650 92 D-01*6.D0,
     5           a33=  0.164 362 784 479 61 D0*6.D0,
     6           a34= -0.135 033 722 413 48 D-01*7.D0,
     7           a35= -0.148 343 453 524 72 D-01*8.D0,
     8           a36=  0.579 229 536 280 84 D-03*9.D0,
     9           a37=  0.323 089 047 037 11 D-02*9.D0,
     1           a38=  0.809 648 029 962 15 D-04*10.D0,
     2           a39= -0.165 576 797 950 37 D-03*10.D0,
     3           a40= -0.449 238 990 618 15 D-04*11.D0)
      PARAMETER (b01=  0.106 580 700 285 13 D01,
     1           b02= -0.157 328 452 902 39 D02,
     2           b03=  0.209 443 969 743 07 D02,
     3           b04= -0.768 677 078 787 16 D01,
     4           b05=  0.261 859 477 879 54 D01,
     5           b06= -0.280 807 811 486 20 D01,
     6           b07=  0.120 533 696 965 17 D01,
     7           b08= -0.845 668 128 125 02 D-02,
     8           b09= -0.126 543 154 777 14 D01,
     9           b10= -0.115 244 078 066 81 D01,
     1           b11=  0.885 210 439 843 18 D0,
     2           b12= -0.642 077 651 816 07 D0,
     3           b13=  0.384 934 601 866 71 D0,
     4           b14= -0.852 147 088 242 06 D0,
     5           b15=  0.489 722 815 418 77 D01,
     6           b16= -0.305 026 172 569 65 D01,
     7           b17=  0.394 205 368 791 54 D-01,
     8           b18=  0.125 584 084 243 08 D0,
     9           b19= -0.279 993 296 987 10 D0,
     1           b20=  0.138 997 995 694 60 D01,
     2           b21= -0.201 899 150 235 70 D01,
     3           b22= -0.821 476 371 739 63 D-02,
     4           b23= -0.475 960 357 349 23 D0,
     5           b24=  0.439 840 744 735 00 D-01,
     6           b25= -0.444 764 354 287 39 D0,
     7           b26=  0.905 720 707 197 33 D0,
     8           b27=  0.705 224 500 879 67 D0,
     9           b28=  0.107 705 126 263 32 D0,
     1           b29= -0.329 136 232 589 54 D0,
     2           b30= -0.508 710 620 411 58 D0,
     3           b31= -0.221 754 008 730 96 D-01,
     4           b32=  0.942 607 516 650 92 D-01,
     5           b33=  0.164 362 784 479 61 D0,
     6           b34= -0.135 033 722 413 48 D-01,
     7           b35= -0.148 343 453 524 72 D-01,
     8           b36=  0.579 229 536 280 84 D-03,
     9           b37=  0.323 089 047 037 11 D-02,
     1           b38=  0.809 648 029 962 15 D-04,
     2           b39= -0.165 576 797 950 37 D-03,
     3           b40= -0.449 238 990 618 15 D-04)
C
      PARAMETER ( tn= 647.096D0,
     2            rhonq= 1.D0/322.0D0,
     3            r= 0.461526D0)
C
      rho= 1.d0/v
      tau= tn/t 
      del= rho*rhonq
C
      tau2= tau*tau
      tau4= tau2*tau2
      tau12= tau4*tau4*tau4
C
      del2= del*del
C
      phd3= a01/del + tau2*(a09 + tau4*a10 + tau12*tau*(a11 +
     1      tau2*(a12 + del*tau4*tau*(a17 + tau4*(a18 + del*(a23 +
     2      del*(a27 + del*(a30 + del*(a33 + del2*(a35 +
     3      del*(a37 + del2*a40))) ))) ))) )) +
     4      del* (a13 + tau2*(a14 + tau4*(a15 + tau*a16)) +
     6      del* (a19 + tau2*(a20 + tau2*(a21 + tau12*a22)) +
     7      del* (a24 + tau2*(a25 + tau2*a26) +
     8      del* (tau*(a28 + tau2*a29) +
     1      del* (a31 + tau2*a32 +
     2      del* (tau2*a34 +
     3      del2*(tau2*a36 +
     4      del* (a38 + tau*a39))) ))) ))
C
      ph3= b01*dlog(del) + b02 + tau*(b03 + tau*(b04 +
     1     tau4*(tau*b05 + tau4*(b06 + tau2*b07 + tau12*(tau*b08 +
     2     del2*(b17 + tau4*(b18 + del*(b23 + del*(b27 +
     3     del*(b30 + del*(b33 + del2*(b35 + del*(b37 +
     4     del2*b40))) ))) ))) ))) ) +
     2     del* (tau2*(b09 + tau4*b10 + tau12*tau*(b11 + tau2*b12)) +
     4     del* (b13 + tau2*(b14 + tau4*(b15 + tau*b16)) +
     6     del* (b19 + tau2*(b20 + tau2*(b21 + tau12*b22)) +
     7     del* (b24 + tau2*(b25 + tau2*b26) +
     8     del* (tau*(b28 + tau2*b29) +
     1     del* (b31 + tau2*b32 +
     2     del* (tau2*b34 +
     3     del2*(tau2*b36 +
     4     del* (b38 + tau*b39))) ))) )))
C
      gvt3n= r*t*(ph3 + del*phd3)
C
      END
C
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is specific volume in m^3/kg.
C  Second argument is temperature in K.
C  Returns specific internal energy in kJ/kg
C
      FUNCTION UVT3N(V,T)
C
      IMPLICIT double precision (A-H,O-Z)
C
      PARAMETER (a03=  0.209 443 969 743 07 D02,
     3           a04= -0.768 677 078 787 16 D01*2.D0,
     4           a05=  0.261 859 477 879 54 D01*7.D0,
     5           a06= -0.280 807 811 486 20 D01*10.D0,
     6           a07=  0.120 533 696 965 17 D01*12.D0,
     7           a08= -0.845 668 128 125 02 D-02*23.D0,
     8           a09= -0.126 543 154 777 14 D01*2.D0,
     9           a10= -0.115 244 078 066 81 D01*6.D0,
     1           a11=  0.885 210 439 843 18 D0*15.D0,
     2           a12= -0.642 077 651 816 07 D0*17.D0,
     4           a14= -0.852 147 088 242 06 D0*2.D0,
     5           a15=  0.489 722 815 418 77 D01*6.D0,
     6           a16= -0.305 026 172 569 65 D01*7.D0,
     7           a17=  0.394 205 368 791 54 D-01*22.D0,
     8           a18=  0.125 584 084 243 08 D0*26.D0,
     1           a20=  0.138 997 995 694 60 D01*2.D0,
     2           a21= -0.201 899 150 235 70 D01*4.D0,
     3           a22= -0.821 476 371 739 63 D-02*16.D0,
     4           a23= -0.475 960 357 349 23 D0*26.D0,
     6           a25= -0.444 764 354 287 39 D0*2.D0,
     7           a26=  0.905 720 707 197 33 D0*4.D0,
     8           a27=  0.705 224 500 879 67 D0*26.D0,
     9           a28=  0.107 705 126 263 32 D0,
     1           a29= -0.329 136 232 589 54 D0*3.D0,
     2           a30= -0.508 710 620 411 58 D0*26.D0,
     4           a32=  0.942 607 516 650 92 D-01*2.D0,
     5           a33=  0.164 362 784 479 61 D0*26.D0,
     6           a34= -0.135 033 722 413 48 D-01*2.D0,
     7           a35= -0.148 343 453 524 72 D-01*26.D0,
     8           a36=  0.579 229 536 280 84 D-03*2.D0,
     9           a37=  0.323 089 047 037 11 D-02*26.D0,
     2           a39= -0.165 576 797 950 37 D-03,
     3           a40= -0.449 238 990 618 15 D-04*26.D0)
C
      PARAMETER ( tn= 647.096D0,
     2            rhonq= 1.D0/322.0D0,
     3            r= 0.461526D0)
C
      rho= 1.d0/v
      tau= tn/t 
      del= rho*rhonq
C
      tau2= tau*tau
      tau4= tau2*tau2
      tau12= tau4*tau4*tau4
C
      del2= del*del
C
      uvt3n= (a03 + tau*(a04 + tau4*(tau*a05 + tau4*(a06 + tau2*a07 +
     1      tau12*(tau*a08 + del2*(a17 + tau4*(a18 + del*(a23 +
     2      del*(a27 + del*(a30 + del*(a33 + del2*(a35 +
     4      del*(a37 + del2*a40))) ))) ))) )))  +
     2      del* (tau*(a09 + tau4*a10 + tau12*tau*(a11 + tau2*a12)) +
     4      del* (tau*(a14 + tau4*(a15 + tau*a16)) +
     6      del* (tau*(a20 + tau2*(a21 + tau12*a22)) +
     7      del* (tau*(a25 + tau2*a26) +
     8      del* (a28 + tau2*a29 +
     1      del* (tau*a32 +
     2      del* (tau*a34 +
     3      del2*(tau*a36 +
     4      del*a39))) ))) )))*tau*t*r
C
      END
C
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is specific volume in m^3/kg.
C  Second argument is temperature in K.
C  Returns isochoric heat capacity in kJ/(kg * K)
C
      FUNCTION CVVT3N(V,T)
C
      IMPLICIT double precision (A-H,O-Z)
C
      PARAMETER (c04= -0.768 677 078 787 16 D01*2.D0,
     4           c05=  0.261 859 477 879 54 D01*42.D0,
     5           c06= -0.280 807 811 486 20 D01*90.D0,
     6           c07=  0.120 533 696 965 17 D01*132.D0,
     7           c08= -0.845 668 128 125 02 D-02*506.D0,
     8           c09= -0.126 543 154 777 14 D01*2.D0,
     9           c10= -0.115 244 078 066 81 D01*30.D0,
     1           c11=  0.885 210 439 843 18 D0*210.D0,
     2           c12= -0.642 077 651 816 07 D0*272.D0,
     4           c14= -0.852 147 088 242 06 D0*2.D0,
     5           c15=  0.489 722 815 418 77 D01*30.D0,
     6           c16= -0.305 026 172 569 65 D01*42.D0,
     7           c17=  0.394 205 368 791 54 D-01*462.D0,
     8           c18=  0.125 584 084 243 08 D0*650.D0,
     1           c20=  0.138 997 995 694 60 D01*2.D0,
     2           c21= -0.201 899 150 235 70 D01*12.D0,
     3           c22= -0.821 476 371 739 63 D-02*240.D0,
     4           c23= -0.475 960 357 349 23 D0*650.D0,
     6           c25= -0.444 764 354 287 39 D0*2.D0,
     7           c26=  0.905 720 707 197 33 D0*12.D0,
     8           c27=  0.705 224 500 879 67 D0*650.D0,
     1           c29= -0.329 136 232 589 54 D0*6.D0,
     2           c30= -0.508 710 620 411 58 D0*650.D0,
     4           c32=  0.942 607 516 650 92 D-01*2.D0,
     5           c33=  0.164 362 784 479 61 D0*650.D0,
     6           c34= -0.135 033 722 413 48 D-01*2.D0,
     7           c35= -0.148 343 453 524 72 D-01*650.D0,
     8           c36=  0.579 229 536 280 84 D-03*2.D0,
     9           c37=  0.323 089 047 037 11 D-02*650.D0,
     3           c40= -0.449 238 990 618 15 D-04*650.D0)
C
      PARAMETER ( tn= 647.096D0,
     2            rhonq= 1.D0/322.0D0,
     3            r= 0.461526D0) 
C
      rho= 1.d0/v
      tau= tn/t 
      del= rho*rhonq
C
      tau2= tau*tau
      tau4= tau2*tau2
      tau12= tau4*tau4*tau4
C
      del2= del*del
C
      phtt3= c04 + tau4*(tau*c05 + tau4*(c06 + tau2*c07 +
     1       tau12*(tau*c08 + del2*(c17 + tau4*(c18 + del*(c23 +
     2       del*(c27 + del*(c30 + del*(c33 + del2*(c35 +
     4       del*(c37 + del2*c40))) ))) ))) )) +
     2       del* (c09 + tau4*c10 + tau12*tau*(c11 + tau2*c12) +
     4       del* (c14 + tau4*(c15 + tau*c16) +
     6       del* (c20 + tau2*(c21 + tau12*c22) +
     7       del* (c25 + tau2*c26 +
     8       del* (tau*c29 +
     1       del* (c32 +
     2       del* (c34 +
     3       del2*c36))) ))) )
C
      cvvt3n= -r*tau*tau*phtt3
C
      END
C
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is pressure in MPa.
C  Second argument is temperature in K.
C  Returns specific enthalpy in kJ/kg
C
C*********************************************************************
      FUNCTION HPT5N(P,T)
C*********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      PARAMETER(A1= .31633380476174D+04, A3= -.30425295300000D+03,
     &          A4= .14381758550000D+04, A5= -.34062035666667D+03,
     &          A6= .34344663500000D+02, A7= -.63611419075956D+01,
     &          A8= .10049580368421D+01, A10=-.16500637020131D-01,
     &          A11=.17887679267013D-03)
C
      TAU=T*1.D-3
      TAUM1=1.D0/TAU
      TAUM2=TAUM1*TAUM1
C
      HPT5N=(((A10*TAUM2*TAUM2*TAUM2+A11*P)*P+A7)*P*TAUM1+A3)*TAUM1+
     &      ((A6*TAU+A5)*TAU+A4)*TAU*TAU+A8*P+A1
C
      END
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is pressure in MPa.
C  Second argument is temperature in K.
C  Returns specific entropy in kJ/(kg*K)
C
C*********************************************************************
      DOUBLE PRECISION FUNCTION SPT5N(P,T)
C*********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
      PARAMETER (A2= -.12459122470994D+02, A3=  .15212647650000D+00,
     &           A4= -.28763517100000D+01, A5=  .51093053500000D+00,
     &           A6= -.45792884666667D-01, A7=  .42407612717304D-02,
     &           A9= -.57982358693700D-04, A10= .14667232906783D-04,
     &           A11=-.11925119511342D-06, R=461.526D-3)

C
      TAU=T*1.D-3
      TAUM1=1.D0/TAU
      TAUM2=TAUM1*TAUM1
      TAUM3=TAUM2*TAUM1
      TAUM6=TAUM3*TAUM3
C     
      SPT5N=-(((A11*P+A10*TAUM6)*P+A7)*P*TAUM3+A9*P+
     &      ((A6*TAU+A5)*TAU+A4)*TAU+A2+A3*TAUM2)-R*DLOG(P*1.D6)
C       
      END
C*********************************************************************
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is pressure in MPa.
C  Second argument is temperature in K.
C  Returns specific volume in m^3/kg
C
C*********************************************************************
      FUNCTION VPT5N(P,T)
C*********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      PARAMETER (A7= -.21203806358652D-02, A8=   .10049580368421D-02,
     &           A9= -.57982358693700D-04, A10= -.36668082266957D-05,
     &           A11= .17887679267013D-06, R1=461.526D-6)
C
      TAU=T*1.D-3
      TAUM2=1.D0/(TAU*TAU)
C
      VPT5N=((A10*TAUM2*TAUM2*TAUM2+A11*P)*P+A7)*TAUM2
     &      +R1*T/P+A9*TAU+A8
C
      END
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is pressure in MPa.
C  Second argument is temperature in K.
C  Returns speed of sound in m/s
C
C*********************************************************************
      DOUBLE PRECISION FUNCTION WPT5N(P,T)
C*********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      parameter(avp10=-.36668082266957D-11,avp11=.35775358534025D-12,
     &       avt7=  .42407612717304D-05,avt9=-.57982358693700D-07,
     &       avt10=.29334465813565D-07,avt11=-.35775358534025D-09,
     &       r=461.526d0,
     &       Av7= -.21203806358652D-02, Av8=   .10049580368421D-02,
     &       Av9= -.57982358693700D-04, Av10= -.36668082266957D-05,
     &       Av11= .17887679267013D-06,
     &       R1=461.526D-6,
     &       A3= .30425295300000D+00, A4=  .28763517100000D+01,
     &       A5=-.10218610700000D+01, A6=  .13737865400000D+00,
     &       A7= .12722283815191D-01, A10= .13200509616104D-03,
     &       A11= -.35775358534025D-06)
C
      TAU=T*1.D-3
      TAUM1=1.D0/TAU
      TAUM2=TAUM1*TAUM1
      TAUM3=TAUM2*TAUM1
      TAUM6=TAUM3*TAUM3
C
      V5=((Av11*P+Av10*TAUM6)*P+Av7)*TAUM2+Av9*TAU+Av8+R1*T/P
C
      CP5=((A11*P+A10*TAUM6)*P+A7)*P*TAUM3+((A6*TAU+A5)*TAU+A4)*TAU
     &        +A3*TAUM2
C
      dvdt=((avt11*p+avt10*taum6)*p+avt7)*taum3+avt9
      dvdt=dvdt+r/p/1.d6
C
      dvdp=(avp11*p+avp10*taum6)*taum2-r*t/p/p/1.d12
C
      CV5=Cp5*1.d3+T*DVDT*DVDT/DVDP
C
      RK=-CP5/CV5*V5/P/DVDP/1.d3
C
      WPT5N=DSQRT(RK*P*1.d6*V5)
C
      END
C*********************************************************************
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is pressure in MPa.
C  Second argument is temperature in K.
C  Returns specific isobaric heat capacity in kJ/(kg*K)
C
C*********************************************************************
      FUNCTION CPPT5N(P,T)
C*********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      PARAMETER(A3= .30425295300000D+00, A4=  .28763517100000D+01,
     &          A5=-.10218610700000D+01, A6=  .13737865400000D+00,
     &          A7= .12722283815191D-01, A10= .13200509616104D-03,
     &          A11= -.35775358534025D-06)
C
      TAU=T*1.D-3
      TAUM1=1.D0/TAU
      TAUM2=TAUM1*TAUM1
C
      CPPT5N=(((A10*TAUM2*TAUM2*TAUM2+A11*P)*P+A7)*P*TAUM1+A3)*TAUM2+
     &       ((A6*TAU+A5)*TAU+A4)*TAU
C
      END
c
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is pressure in MPa.
C  Second argument is temperature in K.
C  Returns specific internal energy in kJ/kg
C
      FUNCTION UPT5N(P,T)
C
      IMPLICIT double precision (A-H,O-Z)
C
      PARAMETER (a02=  0.685 408 416 344 34 D01,
     2           a03=  0.248 051 489 334 66 D-01*3.D0, 
     3           a04= -0.369 015 349 803 33 D0*2.D0,
     4           a05=  0.311 613 182 139 25 D01,
     5           a06= -0.329 616 265 389 17 D0*2.D0)
      PARAMETER (b02=  0.217 746 787 145 71 D-02,
     2           b03= -0.459 428 208 999 10 D-02*3.D0,
     3           b04= -0.397 248 283 595 69 D-05*9.D0,
     4           b05=  0.129 192 282 897 84 D-06*3.D0)
      PARAMETER (c01= -0.125 631 835 895 92 D-03,
     1           c02=  0.217 746 787 145 71 D-02,
     2           c03= -0.459 428 208 999 10 D-02,
     3           c04= -0.397 248 283 595 69 D-05*2.D0,
     4           c05=  0.129 192 282 897 84 D-06*3.D0)
C
      PARAMETER (tn= 1000.0D0,
     1           tnq= 0.001D0,
     2           r= 0.461526D0)
C
      tau= tn/t
      tau2= tau*tau
      tauinv= t*tnq
C
      pi=p
C
      gat5= a02 + tau*a06 + tauinv* tauinv*(a05 + tauinv*(a04 +
     1      tauinv*a03)) +
     2      pi*(b02 + tau2*(b03 + pi*(pi*b05 +
     3      tau2*tau2*tau2*b04)))
C
      gap5= 1.D0/pi + c01 + tau*(c02 + tau2*(c03 + pi*(pi*c05 +
     3      tau2*tau2*tau2*c04)))
C
      upt5n= r*t*(tau*gat5 - (pi*gap5))
C
      END
C
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is pressure in MPa.
C  Second argument is temperature in K.
C  Returns specific Gibbs free energy in kJ/kg.
C
C*********************************************************************
      DOUBLE PRECISION FUNCTION GPT5N(P,T)
C*********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)

      PARAMETER( A1=  .31633380476174D+04, A2= -.12459122470994D+05,
     &             A3= -.15212647650000D+03, A4= -.14381758550000D+04,
     &             A5=  .17031017833333D+03, A6= -.11448221166667D+02,
     &             A7= -.21203806358652D+01, A8=  .10049580368421D+01,
     &             A9= -.57982358693700D-01, A10= -.18334041133478D-2,
     &             A11= .59625597556709D-04, R=461.526D-3)
C
      TAU=T/1.D3
      TAUM1=1.D0/TAU
      TAUM2=TAUM1*TAUM1
      TAUM6=TAUM2*TAUM2*TAUM2
C
      GPT5N=((A11*P+A10*TAUM6)*P+A7)*P*TAUM2+(A9*TAU+A8)*P+
     &      (((A6*TAU+A5)*TAU+A4)*TAU+A2)*TAU+A1+A3*TAUM1
     &      + R*T*DLOG(P*1.D6)       
C       
      END
C*********************************************************************
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is pressure in MPa.
C  Second argument is temperature in K.
C  Returns isochoric heat capacity in kJ/(kg * K)
C
C*********************************************************************
      DOUBLE PRECISION FUNCTION CVPT5N(P,T)
C*********************************************************************
C
      IMPLICIT double precision(A-H,O-Z)
C
      parameter(avp10=-.36668082266957D-11,avp11=.35775358534025D-12,
     &       avt7=  .42407612717304D-05,avt9=-.57982358693700D-07,
     &       avt10=.29334465813565D-07,avt11=-.35775358534025D-09,
     &       r=461.526d0,
     &       A3= .30425295300000D+00, A4=  .28763517100000D+01,
     &       A5=-.10218610700000D+01, A6=  .13737865400000D+00,
     &       A7= .12722283815191D-01, A10= .13200509616104D-03,
     &       A11= -.35775358534025D-06)
C
      TAU=T*1.D-3
      TAUM1=1.D0/TAU
      TAUM2=TAUM1*TAUM1
      TAUM3=TAUM2*TAUM1
      TAUM6=TAUM3*TAUM3
C
      CP5=((A11*P+A10*TAUM6)*P+A7)*P*TAUM3+((A6*TAU+A5)*TAU+A4)*TAU
     &        +A3*TAUM2
C
      dvdt=((avt11*p+avt10*taum6)*p+avt7)*taum3+avt9
      dvdt=dvdt+r/p/1.d6
C
      dvdp=(avp11*p+avp10*taum6)*taum2-r*t/p/p/1.d12
C
      CVPT5N=Cp5+T*DVDT*DVDT/DVDP/1.d3
C
      END
C*********************************************************************
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  Argument is temperature in K.
C  Returns vapour pressure in MPa
C
        function psattn(ts)
        implicit double precision (a-h,o-z)
        parameter(a1=  1167.0521452767d0)
        parameter(a2= -724213.16703206d0)
        parameter(a3= -17.073846940092d0)
        parameter(a4=  12020.824702470d0)
        parameter(a5= -3232555.0322333d0)
        parameter(a6=  14.915108613530d0*2.d0)
        parameter(a7= -4823.2657361591d0*2.d0)
        parameter(a8=  405113.40542057d0*2.d0)
        parameter(a9= -.23855557567849d0)
        parameter(a0=  650.17534844798d0)
        y=ts+a9/(ts-a0)
        b=a5+y*(a4+y*a3)
        c=a8+y*(a7+y*a6)
        ps=c/(dsqrt(b**2-((a1+y)*y+a2)*c*2.d0)-b)
        ps=ps**2
        psattn=ps**2
        end
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is pressure in MPa.
C  Second argument is specific enthalpy in kJ/kg.
C  Returns temperature in K
C
      double precision FUNCTION tph1n(p,h)
C
      IMPLICIT double precision (A-H,O-Z)
C
      PARAMETER (a1 = -.23872489924521D+03,
     1           a2 =  .40421188637945D+03,
     2           a3 =  .11349746881718D+03,
     3           a4 = -.58457616048039D+01,
     4           a5 = -.15285482413140D-03,
     5           a6 = -.10866707695377D-05,
     6           a7 = -.13391744872602D+02,
     7           a8 =  .43211039183559D+02,
     8           a9 = -.54010067170506D+02,
     9           a10=  .30535892203916D+02,
     1           a11= -.65964749423638D+01,
     2           a12=  .93965400878363D-02,
     3           a13=  .11573647505340D-06,
     4           a14= -.25858641282073D-04,
     5           a15= -.40644363084799D-08,
     6           a16=  .66456186191635D-07,
     7           a17=  .80670734103027D-10,
     8           a18= -.93477771213947D-12,
     9           a19=  .58265442020601D-14,
     1           a20= -.15020185953503D-16)
C
      hn = 1.0d0+h*4.d-04
      h2 = hn*hn
      h4 = h2*h2
      h6 = h4*h2
      h10= h6*h4
C
      tph1n =
     1 a1             + hn*(a2          + hn*(a3          +
     2 h4*(a4         + h10*h6*(a5      + h10*(a6         +
     3 p*(a13         + p*(a15          + p*(a17          +
     4 p*(a18         + p*(a19          + p*a20)))))))))) +
     5 p*(a7          + hn*(a8          + hn*(a9          +
     6 hn*(a10        + hn*(a11         + h6*(a12         +
     7 p*(a14         + p*a16)))))))
C
      END
C
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is pressure in MPa.
C  Second argument is specific entropy in kJ/(kg*K).
C  Returns temperature in K
C
      FUNCTION TPS1N(P,S)
C
      IMPLICIT double precision (A-H,O-Z)
C
      PARAMETER (
     &  a1=0.17478268058307D+03,
     &  a2=0.34806930892873D+02,
     &  a3=0.65292584978455D+01,
     &  a4=0.33039981775489D+00,
     &  a5=-.19281382923196D-06,
     &  a6=-.24909197244573D-22,
     &  a7=-.26107636489332D+00,
     &  a8=0.22592965981586D+00,
     &  a9=-.64256463395226D-01,
     &  a10=0.78876289270526D-02,
     &  a11=0.35672110607366D-09,
     &  a12=0.17332496994895D-23,
     &  a13=0.56608900654837D-03,
     &  a14=-.32635483139717D-03,
     &  a15=0.44778286690632D-04,
     &  a16=-.51322156908507D-09,
     &  a17=-.42522657042207D-25,
     &  a18=0.26400441360689D-12,
     &  a19=0.78124600459723D-28,
     &  a20=-.30732199903668D-30)
C
      sn=2.d0+s
      sn2=sn*sn
      sn4=sn2*sn2
      sn8=sn4*sn4
      sn12=sn8*sn4
C
      tps1n = ((((a18*p*sn + a16)*sn8 + a15*sn + a14)*sn + a13)*p +
     &       ((a10*sn + a9)*sn + a8)*sn + a11*sn12 +a7)*p +
     &   ((((((((a20*p + a19)*p*sn + a17)*p + a12)*p + a6)*sn12*sn8 +
     &         a5)*sn8 + a4)*sn + a3)*sn + a2)*sn + a1
C
      END
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is pressure in MPa.
C  Second argument is specific enthalpy in kJ/kg.
C  Returns temperature in K
C
        function tph2n(p,h)
c
        implicit double precision (a-h,o-z)
c
        parameter(a01=1089.8952318288d0)
        parameter(a02=849.51654495535d0)
        parameter(a03=-107.81748091826d0)
        parameter(a04=33.153654801263d0)
        parameter(a05=-7.4232016790248d0)
        parameter(a06=11.765048724356d0)
        parameter(a07=1.8445749355790d0)
        parameter(a08=-4.1792700549624d0)
        parameter(a09=6.2478196935812d0)
        parameter(a10=-17.344563108114d0)
        parameter(a11=-200.58176862096d0)
        parameter(a12=271.96065473796d0)
        parameter(a13=-455.11318285818d0)
        parameter(a14=3091.9688604755d0)
        parameter(a15=252266.40357872d0)
        parameter(a16=-0.61707422868339d-02)
        parameter(a17=-0.31078046629583d0)
        parameter(a18=11.670873077107d0)
        parameter(a19=128127984.04046d0)
        parameter(a20=-985549096.23276d0)
        parameter(a21=2822454697.3002d0)
        parameter(a22=-3594897141.0703d0)
        parameter(a23=1722734991.3197d0)
        parameter(a24=-13551.334240775d0)
        parameter(a25=12848734.664650d0)
        parameter(a26=1.3865724283226d0)
        parameter(a27=235988.32556514d0)
        parameter(a28=-13105236.545054d0)
        parameter(a29=7399.9835474766d0)
        parameter(a30=-551966.97030060d0)
        parameter(a31=3715408.5996233d0)
        parameter(a32=19127.729239660d0)
        parameter(a33=-415351.64835634d0)
        parameter(a34=-62.459855192507d0)
c
        parameter(b01=1489.5041079516d0)
        parameter(b02=743.07798314034d0)
        parameter(b03=-97.708318797837d0)
        parameter(b04=2.4742464705674d0)
        parameter(b05=-0.63281320016026d0)
        parameter(b06=1.1385952129658d0)
        parameter(b07=-0.47811863648625d0)
        parameter(b08=0.85208123431544d-02)
        parameter(b09=0.93747147377932d0)
        parameter(b10=3.3593118604916d0)
        parameter(b11=3.3809355601454d0)
        parameter(b12=0.16844539671904d0)
        parameter(b13=0.73875745236695d0)
        parameter(b14=-0.47128737436186d0)
        parameter(b15=0.15020273139707d0)
        parameter(b16=-0.21764114219750d-02)
        parameter(b17=-0.21810755324761d-01)
        parameter(b18=-0.10829784403677d0)
        parameter(b19=-0.46333324635812d-01)
        parameter(b20=0.71280351959551d-04)
        parameter(b21=0.11032831789999d-03)
        parameter(b22=0.18955248387902d-03)
        parameter(b23=0.30891541160537d-02)
        parameter(b24=0.13555504554949d-02)
        parameter(b25=0.28640237477456d-06)
        parameter(b26=-0.10779857357512d-04)
        parameter(b27=-0.76462712454814d-04)
        parameter(b28=0.14052392818316d-04)
        parameter(b29=-0.31083814331434d-04)
        parameter(b30=-0.10302738212103d-05)
        parameter(b31=0.28217281635040d-06)
        parameter(b32=0.12704902271945d-05)
        parameter(b33=0.73803353468292d-07)
        parameter(b34=-0.11030139238909d-07)
        parameter(b35=-0.81456365207833d-13)
        parameter(b36=-0.25180545682962d-10)
        parameter(b37=-0.17565233969407d-17)
        parameter(b38=0.86934156344163d-14)
c
        parameter(c01=-3236839855524.2d0)
        parameter(c02=7326335090218.1d0)
        parameter(c03=358250899454.47d0)
        parameter(c04=-583401318515.90d0)
        parameter(c05=-10783068217.470d0)
        parameter(c06=20825544563.171d0)
        parameter(c07=610747.83564516d0)
        parameter(c08=859777.22535580d0)
        parameter(c09=-25745.723604170d0)
        parameter(c10=31081.088422714d0)
        parameter(c11=1208.2315865936d0)
        parameter(c12=482.19755109255d0)
        parameter(c13=3.7966001272486d0)
        parameter(c14=-10.842984880077d0)
        parameter(c15=-0.45364172676660d-01)
        parameter(c16=0.14559115658698d-12)
        parameter(c17=0.11261597407230d-11)
        parameter(c18=-0.17804982240686d-10)
        parameter(c19=0.12324579690832d-06)
        parameter(c20=-0.11606921130984d-05)
        parameter(c21=0.27846367088554d-04)
        parameter(c22=-0.59270038474176d-03)
        parameter(c23=0.12918582991878d-02)
c
      if (p.gt.4.d0) goto 20
        p1=p
        p2=p1*p1
        p3=p2*p1
        h1=.0005d0*h-2.1d0
        h2=h1*h1
        h4=h2*h2
        h6=h4*h2
        hh=h6*h6
        tph2n=((( ((( ((( ((( ((( (((
     +  ((a33*p2+a28)*p1+a25)*p1+a23)*h2
     +  +a31*p3+a22)*h2
     +  +a21)*h2
     +  +a20)*h2
     +  +a30*p3+a19)*h2
     +  +a32*p2*p2)*p1*h2
     +  +a15*hh+(a29*p1+a27)*p3)*h4
     +  +a34*p3*p3)*h4
     +  +a24*p2)*h6
     +  +a14)*h6
     +  +a26*p3)*h1
     +  +a13)*h2
     +  +a12)*h2
     +  +a18*p1+a11)*h4
     +  +a10)*h1
     +  +a17*p1+a09)*h1
     +  +a08)*h1
     +  +a16*p1+a07)*p1
     +  +((((a06*h1*hh+a05)*h4+a04)*h1+a03)*h1+a02)*h1+a01
      goto 99
C
 20   P585=(0.12809002730136D-03*H-0.67955786399241D0)*H +
     1       0.90584278514723D+03
c
      if (p.gt.p585) goto 30
        p1=p-2.d0
        p2=p1*p1
        p3=p2*p1
        h1=.0005d0*h-2.6d0
        h2=h1*h1
        h4=h2*h2
        h6=h4*h2
        tph2n=((( ((( ((
     +  ((((b38*p3*p1+b33)*p1+b30)*p2+b20)*p1+b16)*p1+b08)*h6*h6
     +  +(((b36*p1+b34)*p2+b29)*p3+b15)*p1+b07)*h4
     +  +(((b32*p1+b28)*p1+b24)*p2+b14)*p1+b06)*h6
     +  +(((b31*p1+b27)*p2+b19)*p1+b13)*p1+b05)*h6
     +  +((b26*p1+b23)*p2+b12)*p1+b04)*h6
     +  +(b18*p1*h2+b11)*p1)*h4
     +  +((((b35*p3+b25)*p1+b22)*p1+b17)*p1+b10)*p1+b03)*h1
     +  +(b37*p3*p3+b21)*p3+b02)*h1+b09*p1+b01
      goto 99
  30    p1=p+25.d0
        p2=p1*p1
        h1=.0005d0*h-1.8d0
        h2=h1*h1
        h4=h2*h2
        rp1=1.d0/p1
        tph2n=
     +   rp1*(c09+h2*c10
     +  +rp1*(c07+h1*c08
     +  +rp1*rp1*rp1*(c05+h2*c06
     +  +rp1*(c03+h2*c04
     +  +rp1*(c01+h4*c02)))))
     +  +    c11+h1*c12
     +  +p1*(h4*(c13+h4*c14)
     +  +p1*(h4*c15
     +  +p2*p2*(c16+h1*(c17+h2*h1*(c18+h4*h2*
     +          (c19+h2*(c20+h4*(c21+h4*(c22+h2*c23))) )))
     +  )))
 99     end
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  First argument is pressure in MPa.
C  Second argument is specific entropy in kJ/(kg*K).
C  Returns temperature in K
C
      FUNCTION TPS2N(P,S)
c
      IMPLICIT double precision (A-H,O-Z)
C
      parameter (a1 = -0.39235983861984D+06,
     &     a2 =  0.51526573827270D+06,
     &     a3 =  0.40482443161048D+05,
     &     a4=  -0.32193790923902D+03,
     &     a5 =  0.96961424218694D+02,
     &     a6 = -0.22867846371773D+02,
     &     a7 = -0.44942914124357D+06,
     &     a8 = -0.50118336020166D+04,
     &     a9 =  0.35684463560015D+00,
     &     a10=   0.44235335848190D+05,
     &     a11=  -0.13673388811708D+05,
     &     a12=   0.42163260207864D+06,
     &     a13=   0.22516925837475D+05,
     &     a14=   0.47442144865646D+03,
     &     a15=  -0.14931130797647D+03,
     &     a16=  -0.19781126320452D+06,
     &     a17=  -0.23554399470760D+05,
     &     a18=  -0.19070616302076D+05,
     &     a19=   0.55375669883164D+05,
     &     a20=   0.38293691437363D+04,
     &     a21=  -0.60391860580567D+03,
     &     a22=   0.19363102620331D+04,
     &     a23=   0.42660643698610D+04,
     &     a24=  -0.59780638872718D+04,
     &     a25=  -0.70401463926862D+03,
     &     a26=   0.33836784107553D+03,
     &     a27=   0.20862786635187D+02,
     &     a28=   0.33834172656196D-01,
     &     a29=  -0.43124428414893D-04,
     &     a30=   0.16653791356412D+03,
     &     a31=  -0.13986292055898D+03,
     &     a32=  -0.78849547999872D+00,
     &     a33=   0.72132411753872D-01,
     &     a34=  -0.59754839398283D-02,
     &     a35=  -0.12141358953904D-04,
     &     a36=   0.23227096733871D-06,
     &     a37=  -0.10538463566194D+02,
     &     a38=   0.20718925496502D+01,
     &     a39=  -0.72193155260427D-01,
     &     a40=   0.20749887081120D-06,
     &     a41=  -0.18340657911379D-01,
     &     a42=   0.29036272348696D-06,
     &     a43=   0.21037527893619D+00,
     &     a44=   0.25681239729999D-03,
     &     a45=  -0.12799002933781D-01,
     &     a46=  -0.82198102652018D-05)
c
      parameter (b1 = 0.31687665083497D+06,
     &     b2 =  0.20864175881858D+02,
     &     b3 = -0.39859399803599D+06,
     &     b4 = -0.21816058518877D+02,
     &     b5 =  0.22369785194242D+06,
     &     b6 = -0.27841703445817D+04,
     &     b7 =  0.99207436071480D+01,
     &     b8 = -0.75197512299157D+05,
     &     b9 =  0.29708605951158D+04,
     &     b10= -0.34406878548526D+01,
     &     b11=  0.38815564249115D+00,
     &     b12=  0.17511295085750D+05,
     &     b13= -0.14237112854449D+04,
     &     b14=  0.10943803364167D+01,
     &     b15=  0.89971619308495D+00,
     &     b16= -0.33759740098958D+04,
     &     b17=  0.47162885818355D+03,
     &     b18= -0.19188241993679D+01,
     &     b19=  0.41078580492196D+00,
     &     b20= -0.33465378172097D+00,
     &     b21=  0.13870034777505D+04,
     &     b22= -0.40663326195838D+03,
     &     b23=  0.41727347159610D+02,
     &     b24=  0.21932549434532D+01,
     &     b25= -0.10320050009077D+01,
     &     b26=  0.35882943516703D+00,
     &     b27=  0.52511453726066D-02,
     &     b28=  0.12838916450705D+02,
     &     b29= -0.28642437219381D+01,
     &     b30=  0.56912683664855D+00,
     &     b31= -0.99962954584931D-01,
     &     b32= -0.32632037778459D-02,
     &     b33=  0.23320922576723D-03,
     &     b34= -0.15334809857450D+00,
     &     b35=  0.29072288239902D-01,
     &     b36=  0.37534702741167D-03,
     &     b37=  0.17296691702411D-02,
     &     b38= -0.38556050844504D-03,
     &     b39= -0.35017712292608D-04,
     &     b40= -0.14566393631492D-04,
     &     b41=  0.56420857267269D-05,
     &     b42=  0.41286150074605D-07,
     &     b43= -0.20684671118824D-07,
     &     b44=  0.16409393674725D-08)
c
      parameter (c1= 0.90968501005365D+03,
     &     c2 =  0.24045667088420D+04,
     &     c3 = -0.59162326387130D+03,
     &     c4 =  0.54145404128074D+03,
     &     c5 = -0.27098308411192D+03,
     &     c6 =  0.97976525097926D+03,
     &     c7 = -0.46966772959435D+03,
     &     c8 =  0.14399274604723D+02,
     &     c9 = -0.19104204230429D+02,
     &     c10=  0.53299167111971D+01,
     &     c11= -0.21252975375934D+02,
     &     c12= -0.31147334413760D+00,
     &     c13=  0.60334840894623D+00,
     &     c14= -0.42764839702509D-01,
     &     c15=  0.58185597255259D-02,
     &     c16= -0.14597008284753D-01,
     &     c17=  0.56631175631027D-02,
     &     c18= -0.76155864584577D-04,
     &     c19=  0.22440342919332D-03,
     &     c20= -0.12561095013413D-04,
     &     c21=  0.63323132660934D-06,
     &     c22= -0.20541989675375D-05,
     &     c23=  0.36405370390082D-07,
     &     c24= -0.29759897789215D-08,
     &     c25=  0.10136618529763D-07,
     &     c26=  0.59925719692351D-11,
     &     c27= -0.20677870105164D-10,
     &     c28= -0.20874278181886D-10,
     &     c29=  0.10162166825089D-09,
     &     c30= -0.16429828281347D-09)
c
      if(p.gt.4.d0) goto 20
c
      sn=s*.5d0-2.d0
      pn=dsqrt(dsqrt(p))
      pninv=1.d0/pn
C
      sn2=sn*sn
      sn3=sn2*sn
      sn4=sn3*sn
      sn8=sn4*sn4
      sn11=sn8*sn3
C
      sninv=1.d0/sn
      sninv2=sninv*sninv
      sninv4=sninv2*sninv2
      sninv5=sninv4*sninv
      sninv6=sninv5*sninv
      sninv7=sninv6*sninv
      sninv14=sninv7*sninv7
c
      tps2n=((((((((((sninv*a1 + a2)*sninv4 + a3)*sninv6 + a4)*sninv2 +
     &      a5)*sninv + a6)*pninv*sninv4 + (sninv4*a7 + a8)*sninv4*
     &      sninv5 + a9)*pninv + (((((sninv5*a10 + a11)*sninv4 + a12)*
     &      sninv + a13)*sninv7 + a14)*sninv + a15)*sninv2)*pninv +
     &   (sninv*a16 + a17)*sninv4*sninv4)*pninv + (((sn*sninv14*a18 +
     &      a19)*sninv4 + a20)*sninv2 + a21)*sninv)*pninv + ((sninv2*
     &      a22 + a23)*sninv14 + a24)*sninv5 + a25)*pninv*sninv6 +
     & ((((((sn11*sn2*a46 + a45)*pn*sn2 + sn*sn11*a44 + a43)*pn +
     &      (sn11*a42 + a41)*sn4)*pn*sn3 + ((sn8*a40 + a39)*sn*sn4 +
     &      a38)*sn4 + a37)*pn + (((((sn2*a36 + a35)*sn4 + a34)*sn4 +
     &      a33)*sn + a32)*sn4 + a31)*sn + a30)*pn + ((sn3*a29 + a28)*
     &      sn4 + a27)*sn4 + sn*a26)*pn
C
      goto 99
c
 20   if (s.lt.5.85d0) goto 30
c
      sn=10.d0-s*0.12733987011333d+01
      pinv=1.d0/p
      pinvsn=pinv*sn
C
      sn2=sn*sn
      sn3=sn2*sn
      sn4=sn3*sn
C
      TPS2N=((((((sn*b44 + b43)*sn + b42)*p + sn*b41 + b40)*p +
     &      (sn2*b39 + b38)*sn + b37)*p + (sn4*b36 + b35)*sn +
     &      b34)*p + ((((sn*b33 + b32)*sn4 + b31)*sn + b30)*sn +
     &      b29)*sn + b28)*p + (((((sn3*b27 + b26)*sn + b25)*sn +
     &      b24)*sn2 + b23)*sn + b22)*sn +
     &      (((((((((pinv*b2 + b4)*pinv + b7)*pinv + sn*b11 +
     &      b10)*pinvsn + b15)*pinvsn + b20)*sn + b19)*sn3 +
     &      b18)*sn4 + b17)*sn + b16 + ((((pinv*b1 + b3)*pinv +
     &      sn*b6 + b5)*pinv + sn*b9 + b8)*pinv + (sn*sn4*b14 +
     &      b13)*sn + b12)*pinv )*pinv + b21
c
      goto 99
c
 30   continue
      sn=2.d0-s*0.34186865406311d0
      sn2=sn*sn
C
      pinv=1.d0/p
C
      TPS2N=((((((((((sn*c30 + c29)*sn + c28)*sn2 + c27)*sn +
     &      c26)*p + sn*c25 + c24)*p + (sn*c23 + c22)*sn + c21)*
     &      p + (sn*sn2*c20 + c19)*sn + c18)*p + (sn2*sn2*c17 + c16)*
     &      sn + c15)*p + (sn*c14 + c13)*sn + c12)*p + ((sn*c11 + c10)*
     &      sn2 + c9)*sn + c8)*p + ((sn*c7 + c6)*sn + c5)*sn + c4 +
     &      ((sn*c2 + c1)*pinv + c3)*pinv
C
 99   END
c
c
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  Argument is pressure in MPa.
C  Returns saturation temperature in K
C
        function tsatpn(ps)
        implicit double precision (a-h,o-z)
        parameter(a1=  1167.0521452767d0)
        parameter(a2= -724213.16703206d0*(-1.d0))
        parameter(a3= -17.073846940092d0)
        parameter(a4=  12020.824702470d0)
        parameter(a5= -3232555.0322333d0*(-1.d0))
        parameter(a6=  14.915108613530d0)
        parameter(a7= -4823.2657361591d0)
        parameter(a8=  405113.40542057d0*(-1.d0))
        parameter(a9= -.23855557567849d0)
        parameter(a0=  650.17534844798d0*.5d0)
        parameter(aa=-4.d0*a0)
        parameter(ab= 4.d0*a0**2-a9)
        y=dsqrt(dsqrt(ps))
        f=a7+y*(a4+y*a1)
        g=a8+y*(a5+y*a2)
        t1=g/(f+dsqrt(f**2+((a3+y)*y+a6)*g*4.d0))+a0
        tsatpn=t1-dsqrt((aa+t1)*t1+ab)
        end
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  The argument is pressure in MPa.
C  Returns temperature in K.
C
      FUNCTION FB23P(P)
C
C  function for the calculation of the temperature at the boundary
C  between the regions 2 and 3
C
      IMPLICIT double precision (A-H,O-Z)
C
      FB23P = 0.57254459862746D+03 + DSQRT((P - 0.13918839778870D+02)
     &        *0.98106832075625D+03)
C
      END
c
c                 IAPWS Industrial Formulation 1997
c        for the Thermodynamic Properties of Water and Steam
c                           (IAPWS-IF97)
C
C  The argument is temperature in K.
C  Returns pressure in MPa.
C
      FUNCTION FB23(T)
C
      IMPLICIT double precision (A-H,O-Z)
C
      FB23=(0.10192970039326D-02*T - 0.11671859879975D+01)*T +
     1      0.34805185628969D+03
C
      END
C
C*********************************************************************
      BLOCKDATA ZUSPKT
C*********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      COMMON/CCRIT/VC,SC,HC
      COMMON/CGR/PGR,TGR13
      DATA VC  /3.105590062111801D-3/
     &     SC  /4.41202148223476D0/
     &     HC  /2087.54684511650D0/
      DATA PGR  /16.52916425260452D0/
     &     TGR13 /623.15D+0/
      END
C*********************************************************************
      BLOCKDATA  SUBST
C*********************************************************************
      IMPLICIT double precision(A-H,O-Z)
      COMMON/CSUB2/TC,PC,DC
      DATA TC   /0.647096000D+03/
     *     PC   /0.220640000D+02/
     *     DC   /0.322000000D+03/
      END
