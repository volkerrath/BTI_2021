%  
% Path ->  /home/volker/Matlab/heat1d_toolbox
% 
%   wm                  - wm(J) calculates the diagonal weighting 
%   get_model           - constructs model parameter index ip
%   paleo_grid          - defines temporal inversion grid 
%   put_props           - writes geological units and associated physical parameter 
%   reg1d               - generates the regularization matrix L 
%   c2n                 - intepolates cell-to-node linearly
%   m2p                 - from inverse parameter space m to physical parameter p
%   rhofT               - calculate the density i(in kg/m^3) of pure water, 
%   rhoiT               - ice density in [kg/m**3]
%   n2c                 - intepolates cell-to-node linearly
%   kfT                 - calculate the thermal conductivity kf in W/(m*K) of
%   kiT                 - ice thermal conductivity [mW/(m*K)]
%   p2m                 - from physical parameter space p to inverse parameter m
%   kmT                 - thermal conductivity as function of temperature
%   rms                 - calculates the rms value
%   tri                 - window.
%   paleo_osterkamp     - initializes paleotemperature of permfrost surface T in
%   paleo_step          - initializes step function for paleoclimate
%   threshsp            - contracts full (covariance) matrix to sparse matrix S
%   cglsACB             - applies conjugate gradient algorithm 
%   cglssh              - Conjugate gradient algorithm applied implicitly to the normal equations.
%   paleo_boxcar        - initializes step function for paleoclimate
%   heat1ds             - solves nonlinear steady-state heat equation
%   heat1dt             - solves nonlinear time-dependent heat equation
%   get_data            - constructs inversion data from temperature logs
%   mcgls               - 
%   set_cell            - (No help available)
%   paleo_boxcar_smooth - initializes step function for paleoclimate
%   cgls                - applies conjugate gradient algorithm 
%   cpfT                - pure water heat capacity depending on
%   cpiT                - ice isobaric heat capacity [J/(kg*K]
%   cpmT                - heat capacity of rocks as function of temperature
%   smooth1             - smoothes input time series f with triangular  weights
%   sensfds_pet         - calculates stationary nonlinear Jacobian 
%   set_mesh            - generates a mesh given type
%   jacfds_pet          - calculates transient Jacobians with respect to paleoclimate. 
%   get_props           - reads geological units and associated physical parameter 
%   heat1dns            - solves nonlinear stationary heat equation
%   heat1dnt            - solves nonlinear time-dependent heat equation 
%   heat1dsg            - solves nonlinear steady-state heat equation
%   heat1dtg            - solves nonlinear time-dependent heat equation
%   l_1d                - generates the regularization matrix L
%   paleo_rellstab      - initializes paleoclimate after Haenel (1988)
%   set_gst_prior       - sets GST prior
%   paleo_haenel        - initializes paleoclimate after Haenel (1988)
%   sensfdt_pal         - calculates transient Jacobians with respect to paleoclimate 
%   paleo_transylvania  - initializes step function for paleoclimate
%   jacfdt_pal          - calculates transient Jacobians with respect to paleoclimate. 
%   jacfdt_pet          - calculates transient Jacobians with respect to petrophysical
%   cov1_pp             - calculates a gaussian covariance a-priori matrix
%   wrobust             - calculates IRLS weights for robust regression 
%   we1d                - we1d(m,x,mode,eps,bet)
%   heat1dans           - solves analytically the stationary heat equation 
%   heat1dant           - solves analytically time-dependent heat equation (N&B,89)
%   put_data            - writes DATA and ERRORS 
%   set_paleo_grid      - defines temporal inversion grid 
%   wmean               - weighted mean
%   ftheta              - calculates the fluid/ice  
%   plot_d              - plots (multiple) data
%   plot_m              - plots parameter + error bars 
%   wp_ms               - the minimum support (MS) weighting function
%   paleo_huang         - initializes global warming after Huang & Pollack (2001) 
