	/*****************************************************************              
	 *  INSTITUT FUER ENERGIE- UND KLIMAFORSCHUNG                    *              
	 +  IEK-5 PHOTOVOLTAIK                                           *              
	 *                                                               *              
	 *        ########                _   _                          *              
	 *     ##########                |_| |_|                         *              
	 *    ##########     ##         _ _   _ _     ___ ____ _   _     *              
	 *   ##########     ####       | | | | | |   |_ _/ ___| | | |    *              
	 *   #########     #####    _  | | | | | |    | | |   | |_| |    *              
	 *   #    ###     ######   | |_| | |_| | |___ | | |___|  _  |    *              
	 *    #          ######     \___/ \___/|_____|___\____|_| |_|    *              
	 *     ##      #######      F o r s c h u n g s z e n t r u m    *              
	 *       ##########                                              *              
	 *                                                               *              
	 *   http://www.fz-juelich.de/iek/iek-5/DE/Home/home_node.html   *              
	 *****************************************************************
	 *                                                               *
	 *    Dr. Bart E. Pieters 2015                                   *
	 *                                                               *             
	 *****************************************************************/              
 
PhotoVoltaic MOdule Simulator (PVMOS) version 0.73
=================================================

PVMOS is an ordinary differential equation solver using finite-differences specifically
designed to electrically model solar modules. A paper about PVMOS has been submitted to 
IEEE Journal of photovoltaics [1]. Untill that paper is out I refer the user to [2], 
which is an older version 0.73
older reference which is still relevant is an paper from 2011 which discusse the same 
method as PVMOS uses (Network Simulation Method with a variable and daptive mesh) [3].
There is a first drafdt of a manual which at the moment not much more than a function 
reference. However, in combination with the provided examples it may be enough to get 
you started.

Note that this is still beta software at its pre 1.0 state (i.e. there is still lots of 
features missing that are planned. This may also mean that e.g. syntax of input files 
changes with a new release, be prepared to rewrite input files if you want to stay up 
to date.).



[1] B. E. Pieters, "PVMOS: A Free and Open Source Simulation Tool for Solar Modules," 
    submitted to Journal of Photovoltaics, 2014.
    
[2] B. E. Pieters, "A free and open source finite-difference simulation tool for solar 
    modules." Photovoltaic Specialist Conference (PVSC), 2014 IEEE 40th. IEEE, 2014.
    
[3] B. E. Pieters, "Spatial modeling of Thin-Film Solar Modules using the Network 
    Simulation Method and SPICE," Journal of Photovoltaics, 2011.
