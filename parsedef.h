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
 *    Dr. Bart E. Pieters 2015                                   *
 *                                                               *             
 *****************************************************************/                                                                             
/* parse keys and keywords */

/* enummeration of keys */
typedef enum {
	NEWMESH,
	ADDCOL_E,
	ADDROW_N,
	ADDCOL_W,
	ADDROW_S,
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
	MOVEMESH,
	ROTATEMESH,
	FLIPX,
	FLIPY,
	SCALE,
	SCALEX,
	SCALEY,
	SETBB,
	RESOLVPOLY,
	RESOLVCIRC,
	RESOLVRECT,
	SIMPLIFY_MESH,
	LOADMESH,
	SAVEMESH,
	PRINTMESH,
	PRINTCONN,
	PRINTSURF,
	PRINTPOT,
	PRINTIV,
	PRINTSOLPAR,
	PRINTPROBE,
	PRINTINIPIV,
	PRINTPARS,
	PRINTLOCALJV,
	SURFCOLCUR,
	SURFDCOLCUR,
	SURFVPLOT,
	SURFPPLOT,
	SURFJPLOT,
	SURFVJPLOT,
	SURFEPLOT,
	LOAD_POLY,
	DEF_POLY,
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
	SET_PTJV,
	SET_PTEXTJV,
	SET_R,
	SET_T,
	SET_SPLITY,			
	SET_SPLITX,		
	SET_SEL_REL,
	SET_SEL_RVP,
	SET_SEL_RVN,
	SET_SEL_JV,
	SET_SEL_2DJV,
	SET_SEL_1DJV,
	SET_SEL_PTJV,
	SET_SEL_PTEXTJV,
	SET_SEL_R,
	SET_SEL_T,
	SET_SEL_SPLITY,			
	SET_SEL_SPLITX,	
	MAXITER,
	TOLV,	
	RELTOLV,
	TOLKCL,	
	RELTOLKCL,
	NLINSEARCH,
	GMINSTEP,
	GMINSTART,
	GMINMAX,
	GMINFAC,
	GMIN,				
	SOLVE,					
	REFINEOC,				
	REFINEMPP,			
	ADAPTIVE_SOLVE,
	TIC,
	TOC,
	EXPR_DEF,
	EXPR_DEF_SOLPAR,
	EXPR_DEF_BB,	
	WHILE,		
	ENDWHILE,	
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
      	{"addcol_right", ADDCOL_E},
      	{"addrow_top", ADDROW_N},
      	{"addcol_left", ADDCOL_W},
      	{"addrow_bottom", ADDROW_S},
      	{"rm", RMMESH},
	{"joinmesh", JOINMESH},
	{"joinmesh_h", JOINMESH_H},
	{"joinmesh_v", JOINMESH_V},
	{"dupmesh", DUPMESH},
	{"add_electrode", ADDEL},
	{"split_x", SPLITX},
	{"split_y", SPLITY},
	{"split_xy", SPLITXY},
	{"split_long", SPLITLONG},
	{"split_coarse", SPLITCOARSE},
	{"move", MOVEMESH},
	{"rotate", ROTATEMESH},
	{"flipx", FLIPX},
	{"flipy", FLIPY},
	{"scale", SCALE},
	{"scalex", SCALEX},
	{"scaley", SCALEY},
	{"boundingbox", SETBB},
	{"resolve_poly", RESOLVPOLY},
	{"resolve_circ", RESOLVCIRC},
	{"resolve_rect", RESOLVRECT},
	{"simplify", SIMPLIFY_MESH},
	{"loadmesh", LOADMESH},
	{"savemesh", SAVEMESH},
	{"printmesh", PRINTMESH},
	{"printconn", PRINTCONN},
	{"printarea", PRINTSURF},
	{"printV", PRINTPOT},
	{"printIV", PRINTIV},
	{"print_solpar", PRINTSOLPAR},
	{"print_probe", PRINTPROBE},
	{"printInIp", PRINTINIPIV},
	{"printpars", PRINTPARS},
	{"printlocalJV", PRINTLOCALJV},
	{"surfFc", SURFCOLCUR},
	{"surffc", SURFDCOLCUR},
	{"surfVplot", SURFVPLOT},
	{"surfPplot", SURFPPLOT},
	{"surfJplot", SURFJPLOT},
	{"surfVjplot", SURFVJPLOT},
	{"surfEplot", SURFEPLOT},
	{"load_poly", LOAD_POLY},
	{"define_poly", DEF_POLY},
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
	{"set_PTJV", SET_PTJV},
	{"set_PTJV_extra", SET_PTEXTJV},
	{"set_R", SET_R},
	{"set_T", SET_T},
	{"set_SplitX", SET_SPLITX},			
	{"set_SplitY", SET_SPLITY},		
	{"set_sel_Rel", SET_SEL_REL},
	{"set_sel_Rvp", SET_SEL_RVP},
	{"set_sel_Rvn", SET_SEL_RVN},
	{"set_sel_JV", SET_SEL_JV},
	{"set_sel_2DJV", SET_SEL_2DJV},
	{"set_sel_1DJV", SET_SEL_1DJV},
	{"set_sel_PTJV", SET_SEL_PTJV},
	{"set_sel_PTJV_extra", SET_SEL_PTEXTJV},
	{"set_sel_R", SET_SEL_R},
	{"set_sel_T", SET_SEL_T},
	{"set_sel_SplitX", SET_SEL_SPLITX},			
	{"set_sel_SplitY", SET_SEL_SPLITY},		
	{"maxiter", MAXITER},
	{"tol_V",TOLV},	
	{"rel_tol_V",RELTOLV},
	{"tol_kcl",TOLKCL},	
	{"rel_tol_kcl",RELTOLKCL},
	{"N_Lin_Search",NLINSEARCH},	
	{"GminStep",GMINSTEP},		
	{"GminStart",GMINSTART},	
	{"GminMax",GMINMAX},
	{"GminFac",GMINFAC},		
	{"Gmin",GMIN},	
	{"solve", SOLVE},				
	{"refine_oc", REFINEOC},				
	{"refine_mpp", REFINEMPP},
	{"adaptive_solve", ADAPTIVE_SOLVE},
	{"tic", TIC},
	{"toc", TOC},
	{"define", EXPR_DEF},
	{"define_solpar", EXPR_DEF_SOLPAR},
	{"define_bb", EXPR_DEF_BB},
	{"while", WHILE},
	{"endwhile", ENDWHILE},
	{"out_quiet",_QUIET},
	{"out_normal",_NORMAL},
	{"out_verbose",_VERBOSE},
	{"out_debug",_DEBUG},
      	{NULL, NONE}
};

