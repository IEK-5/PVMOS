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
	RMMESH,
	JOINMESH,
	JOINMESH_H,
	JOINMESH_V,
	DUPMESH,
	ADDEL,
	SPLITX,
	SPLITY,
	SPLITXY,
	SPLITLONG,
	SPLITCOARSE,
	SIMPLIFY_MESH,
	LOADMESH,
	SAVEMESH,
	PRINTMESH,
	PRINTCONN,
	PRINTSURF,
	PRINTPOT,
	PRINTIV,
	PRINTPARS,
	SURFVPLOT,
	SURFPPLOT,
	LOAD_POLY,
	SELECT_RECT,
	SELECT_RECT_CONTOUR,
	SELECT_CIRC,
	SELECT_CIRC_CONTOUR,
	SELECT_POLY,
	SELECT_POLY_CONTOUR,
	SELECT_AREA,
	DESELECT,
	ASSIGN_PROP,
	SET_REL,
	SET_RVP,
	SET_RVN,
	SET_JV,
	SET_2DJV,
	SET_1DJV,
	SET_R,
	SET_T,
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
      	{"rm", RMMESH},
	{"joinmesh", JOINMESH},
	{"joinmesh_h", JOINMESH_H},
	{"joinmesh_v", JOINMESH_V},
	{"duplicate_mesh", DUPMESH},
	{"add_electrode", ADDEL},
	{"split_x", SPLITX},
	{"split_y", SPLITY},
	{"split_xy", SPLITXY},
	{"split_long", SPLITLONG},
	{"split_coarse", SPLITCOARSE},
	{"simplify", SIMPLIFY_MESH},
	{"loadmesh", LOADMESH},
	{"savemesh", SAVEMESH},
	{"printmesh", PRINTMESH},
	{"printconn", PRINTCONN},
	{"printarea", PRINTSURF},
	{"printV", PRINTPOT},
	{"printIV", PRINTIV},
	{"printpars", PRINTPARS},
	{"surfVplot", SURFVPLOT},
	{"surfPplot", SURFPPLOT},
	{"load_poly", LOAD_POLY},
	{"select_rect", SELECT_RECT},
	{"select_rect_contour", SELECT_RECT_CONTOUR},
	{"select_circ", SELECT_CIRC},
	{"select_circ_contour", SELECT_CIRC_CONTOUR},
	{"select_poly", SELECT_POLY},
	{"select_poly_contour", SELECT_POLY_CONTOUR},
	{"select_area", SELECT_AREA},
	{"deselect", DESELECT},
	{"assign_properties", ASSIGN_PROP},
	{"set_Rel", SET_REL},
	{"set_Rvp", SET_RVP},
	{"set_Rvn", SET_RVN},
	{"set_JV", SET_JV},
	{"set_2DJV", SET_2DJV},
	{"set_1DJV", SET_1DJV},
	{"set_R", SET_R},
	{"set_T", SET_T},
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

