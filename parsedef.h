/*
    PhotoVoltaic Module Simulator (PVMOS), a finite difference ODE solver 
    for solar modules. 
    Copyright (C) 2014  B. E. Pieters, 
    IEK-5 Photovoltaik, Forschunszentrum Juelich

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
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
 *    Dr. Bart E. Pieters 2014                                   *
 *                                                               *             
 *****************************************************************/                                                                             
/* parse keys and keywords */

/* enummeration of keys */
typedef enum {
	NEWMESH,
	JOINMESH,
	JOINMESH_H,
	JOINMESH_V,
	SPLITX,
	SPLITY,
	SPLITXY,
	SPLITLONG,
	SIMPLIFY_MESH,
	LOADMESH,
	SAVEMESH,
	PRINTMESH,
	PRINTCONN,
	PRINTSURF,
	PRINTPOT,
	PRINTMESHSEL,
	PRINTCONNSEL,
	PRINTSURFSEL,
	PRINTPOTSEL,
	PRINTIV,
	PRINTPARS,
	SURFVPLOT,
	SURFPPLOT,
	LOAD_POLY,
	SELECT_RECT,
	SELECT_CIRC,
	SELECT_POLY,
	DESELECT,
	ASSIGN_PROP,
	SET_RP,
	SET_RN,
	SET_RPVP,
	SET_RNVP,
	SET_RPVN,
	SET_RNVN,
	SET_JV,
	SET_2DJV,
	SET_1DJV,
	SET_R,
	SET_SPLITY,			
	SET_SPLITX,		
	MAXITER,
	TOLV,	
	RELTOLV,
	TOLKCL,	
	RELTOLKCL,				
	SOLVE,			
	ADAPTIVE_SOLVE,
	_QUIET,
	_NORMAL,
	_VERBOSE,
	_DEBUG,
	NONE
} PRSDEF;
/* Keyword to key mapping struct */
typedef struct {
	const char *name;
	PRSDEF PAR;
} KeyWord;


/* The keyword table mapping keywords to keys */
const KeyWord KeyTable[] =
{
      	{"newmesh", NEWMESH},
	{"joinmesh", JOINMESH},
	{"joinmesh_h", JOINMESH_H},
	{"joinmesh_v", JOINMESH_V},
	{"split_x", SPLITX},
	{"split_y", SPLITY},
	{"split_xy", SPLITXY},
	{"split_long", SPLITLONG},
	{"simplify", SIMPLIFY_MESH},
	{"loadmesh", LOADMESH},
	{"savemesh", SAVEMESH},
	{"printmesh", PRINTMESH},
	{"printconn", PRINTCONN},
	{"printarea", PRINTSURF},
	{"printV", PRINTPOT},
	{"printmesh_sel", PRINTMESHSEL},
	{"printconn_sel", PRINTCONNSEL},
	{"printarea_sel", PRINTSURFSEL},
	{"printV_sel", PRINTPOTSEL},
	{"printIV", PRINTIV},
	{"printpars", PRINTPARS},
	{"surfVplot", SURFVPLOT},
	{"surfPplot", SURFPPLOT},
	{"load_poly", LOAD_POLY},
	{"select_rect", SELECT_RECT},
	{"select_circ", SELECT_CIRC},
	{"select_poly", SELECT_POLY},
	{"deselect", DESELECT},
	{"assign_properties", ASSIGN_PROP},
	{"set_Rp", SET_RP},
	{"set_Rn", SET_RN},
	{"set_Rpvp", SET_RPVP},
	{"set_Rnvp", SET_RNVP},
	{"set_Rpvn", SET_RPVN},
	{"set_Rnvn", SET_RNVN},
	{"set_JV", SET_JV},
	{"set_2DJV", SET_2DJV},
	{"set_1DJV", SET_1DJV},
	{"set_R", SET_R},
	{"set_SplitX", SET_SPLITX},			
	{"set_SplitY", SET_SPLITY},		
	{"maxiter", MAXITER},
	{"tol_V",TOLV},	
	{"rel_tol_V",RELTOLV},
	{"tol_kcl",TOLKCL},	
	{"rel_tol_kcl",RELTOLKCL},	
	{"solve", SOLVE},
	{"adaptive_solve", ADAPTIVE_SOLVE},
	{"out_quiet",_QUIET},
	{"out_normal",_NORMAL},
	{"out_verbose",_VERBOSE},
	{"out_debug",_DEBUG},
      	{NULL, NONE}
};

