/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__Ih
#define _nrn_initial _nrn_initial__Ih
#define nrn_cur _nrn_cur__Ih
#define _nrn_current _nrn_current__Ih
#define nrn_jacob _nrn_jacob__Ih
#define nrn_state _nrn_state__Ih
#define _net_receive _net_receive__Ih 
#define rates rates__Ih 
#define states states__Ih 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg(int);
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define gbar _p[0]
#define gbar_columnindex 0
#define i _p[1]
#define i_columnindex 1
#define g _p[2]
#define g_columnindex 2
#define tauh _p[3]
#define tauh_columnindex 3
#define ainf _p[4]
#define ainf_columnindex 4
#define a _p[5]
#define a_columnindex 5
#define b _p[6]
#define b_columnindex 6
#define binf _p[7]
#define binf_columnindex 7
#define Da _p[8]
#define Da_columnindex 8
#define Db _p[9]
#define Db_columnindex 9
#define v _p[10]
#define v_columnindex 10
#define _g _p[11]
#define _g_columnindex 11
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_a_tauh(void);
 static void _hoc_a_inf(void);
 static void _hoc_rates(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_Ih", _hoc_setdata,
 "a_tauh_Ih", _hoc_a_tauh,
 "a_inf_Ih", _hoc_a_inf,
 "rates_Ih", _hoc_rates,
 0, 0
};
#define a_tauh a_tauh_Ih
#define a_inf a_inf_Ih
 extern double a_tauh( _threadargsprotocomma_ double );
 extern double a_inf( _threadargsprotocomma_ double );
 /* declare global and static user variables */
#define Vmid_ac Vmid_ac_Ih
 double Vmid_ac = -92;
#define eh eh_Ih
 double eh = -40;
#define h h_Ih
 double h = 0;
#define k_ac k_ac_Ih
 double k_ac = -7.25;
#define m m_Ih
 double m = 1;
#define q10 q10_Ih
 double q10 = 1.5;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "eh_Ih", "mV",
 "Vmid_ac_Ih", "mV",
 "k_ac_Ih", "mV",
 "gbar_Ih", "pS/microm2",
 "i_Ih", "mA/cm2",
 "g_Ih", "pS/microm2",
 "tauh_Ih", "ms",
 "ainf_Ih", "1",
 0,0
};
 static double a0 = 0;
 static double b0 = 0;
 static double delta_t = 0.01;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "eh_Ih", &eh_Ih,
 "Vmid_ac_Ih", &Vmid_ac_Ih,
 "k_ac_Ih", &k_ac_Ih,
 "m_Ih", &m_Ih,
 "h_Ih", &h_Ih,
 "q10_Ih", &q10_Ih,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(NrnThread*, _Memb_list*, int);
static void nrn_state(NrnThread*, _Memb_list*, int);
 static void nrn_cur(NrnThread*, _Memb_list*, int);
static void  nrn_jacob(NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(NrnThread*, _Memb_list*, int);
static void _ode_matsol(NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[0]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"Ih",
 "gbar_Ih",
 0,
 "i_Ih",
 "g_Ih",
 "tauh_Ih",
 "ainf_Ih",
 0,
 "a_Ih",
 "b_Ih",
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 12, _prop);
 	/*initialize range parameters*/
 	gbar = 3;
 	_prop->param = _p;
 	_prop->param_size = 12;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 1, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _IHnew_reg() {
	int _vectorized = 1;
  _initlists();
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 12, 1);
  hoc_register_dparam_semantics(_mechtype, 0, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 Ih /Users/cyprien/Desktop/EBrains_2025/Code/Mechanism/IHnew.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(_threadargsprotocomma_ double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[1], _dlist1[1];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {int _reset = 0; {
   rates ( _threadargscomma_ v ) ;
   Da = ( ainf - a ) / tauh ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
 rates ( _threadargscomma_ v ) ;
 Da = Da  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tauh )) ;
  return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) { {
   rates ( _threadargscomma_ v ) ;
    a = a + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tauh)))*(- ( ( ( ainf ) ) / tauh ) / ( ( ( ( - 1.0 ) ) ) / tauh ) - a) ;
   }
  return 0;
}
 
double a_inf ( _threadargsprotocomma_ double _lV ) {
   double _la_inf;
 _la_inf = 1.0 / ( 1.0 + exp ( - ( _lV - Vmid_ac ) / k_ac ) ) ;
   
return _la_inf;
 }
 
static void _hoc_a_inf(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  a_inf ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double a_tauh ( _threadargsprotocomma_ double _lV ) {
   double _la_tauh;
  _la_tauh = 556.0 + 1100.0 * exp ( - 0.5 * pow( ( ( _lV ) / 11.06 ) , 2.0 ) ) ;
    
return _la_tauh;
 }
 
static void _hoc_a_tauh(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  a_tauh ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
static int  rates ( _threadargsprotocomma_ double _lV ) {
   double _lqt ;
 _lqt = pow( q10 , ( ( celsius - 24.0 ) / 10.0 ) ) ;
   tauh = a_tauh ( _threadargscomma_ _lV ) / _lqt ;
   ainf = a_inf ( _threadargscomma_ _lV ) ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 rates ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 1;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 (_p, _ppvar, _thread, _nt);
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 1; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
  a = a0;
  b = b0;
 {
   rates ( _threadargscomma_ v ) ;
   a = ainf ;
   b = binf ;
   }
 
}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
 initmodel(_p, _ppvar, _thread, _nt);
}
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   g = ( gbar * ( pow( a , m ) ) * ( pow( b , h ) ) ) ;
   i = ( 0.0001 ) * g * ( v - eh ) ;
   }
 _current += i;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
 	}
 _g = (_g - _rhs)/.001;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
 {   states(_p, _ppvar, _thread, _nt);
  }}}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = a_columnindex;  _dlist1[0] = Da_columnindex;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/cyprien/Desktop/EBrains_2025/Code/Mechanism/IHnew.mod";
static const char* nmodl_file_text = 
  "COMMENT \n"
  "\n"
  "Model for an Ih amendola\n"
  "ENDCOMMENT\n"
  "\n"
  "\n"
  "NEURON\n"
  " {\n"
  "  SUFFIX Ih\n"
  "NONSPECIFIC_CURRENT i\n"
  " \n"
  "RANGE gbar, g, i  \n"
  "\n"
  "RANGE tauh\n"
  "RANGE ainf\n"
  "\n"
  "\n"
  " }\n"
  "\n"
  "UNITS {\n"
  "\n"
  "(pS) =(picosiemens)\n"
  "(mV) = (millivolt)\n"
  "(mA) = (milliamp)\n"
  "\n"
  "}\n"
  "\n"
  "PARAMETER { \n"
  "  gbar = 3(pS/microm2)\n"
  "  eh = -40 (mV)\n"
  "  Vmid_ac = -92 (mV)\n"
  "  k_ac = -7.25 (mV) \n"
  "  m=1         \n"
  "  h=0        : gate parameters according to the HH formalism (m*m*m*h)\n"
  "celsius = 32 (degC)\n"
  "q10=1.5\n"
  "}\n"
  "\n"
  " ASSIGNED {\n"
  "  v	(mV)\n"
  "   i 	(mA/cm2)\n"
  "  g	(pS/microm2)\n"
  " tauh (ms)\n"
  "   ainf (1)\n"
  "binf(1)\n"
  " \n"
  " }\n"
  "\n"
  "\n"
  "STATE {a b}\n"
  "\n"
  "BREAKPOINT {\n"
  "  SOLVE states METHOD cnexp\n"
  "  g = (gbar*(a^m)*(b^h)) \n"
  "  i = (0.0001)*g*(v-eh)\n"
  "  \n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "\n"
  "rates(v)\n"
  "a= ainf\n"
  "b=binf\n"
  "\n"
  "\n"
  "}\n"
  "\n"
  "DERIVATIVE states {\n"
  " rates(v)\n"
  "  a' = (ainf-a)/tauh\n"
  " \n"
  "  \n"
  "}\n"
  "\n"
  "\n"
  "FUNCTION a_inf (V (mV)) () {\n"
  "\n"
  "  a_inf = 1/(1+exp(-(V-Vmid_ac)/k_ac))  : activation system \n"
  "}\n"
  "\n"
  "\n"
  "\n"
  "FUNCTION a_tauh (V (mV)) (ms) {\n"
  "UNITSOFF\n"
  "\n"
  "a_tauh= 556+ 1100*exp(-0.5*((V)/11.06)^2)\n"
  ": time constant of activation depends on V \n"
  "\n"
  "UNITSON\n"
  "}\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "PROCEDURE rates(V (mV)) {\n"
  "\n"
  "LOCAL qt\n"
  "\n"
  "qt=q10^((celsius-24)/10)\n"
  "\n"
  "tauh=a_tauh(V)/qt\n"
  "\n"
  "ainf=a_inf(V)\n"
  "\n"
  "\n"
  "}\n"
  ;
#endif
