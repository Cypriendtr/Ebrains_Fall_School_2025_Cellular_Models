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
 
#define nrn_init _nrn_init__kaDasoma
#define _nrn_initial _nrn_initial__kaDasoma
#define nrn_cur _nrn_cur__kaDasoma
#define _nrn_current _nrn_current__kaDasoma
#define nrn_jacob _nrn_jacob__kaDasoma
#define nrn_state _nrn_state__kaDasoma
#define _net_receive _net_receive__kaDasoma 
#define rates rates__kaDasoma 
#define states states__kaDasoma 
 
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
#define taurecov _p[1]
#define taurecov_columnindex 1
#define i _p[2]
#define i_columnindex 2
#define g _p[3]
#define g_columnindex 3
#define atau _p[4]
#define atau_columnindex 4
#define btau _p[5]
#define btau_columnindex 5
#define ainf _p[6]
#define ainf_columnindex 6
#define binf _p[7]
#define binf_columnindex 7
#define a _p[8]
#define a_columnindex 8
#define b _p[9]
#define b_columnindex 9
#define b2 _p[10]
#define b2_columnindex 10
#define ik _p[11]
#define ik_columnindex 11
#define ek _p[12]
#define ek_columnindex 12
#define Da _p[13]
#define Da_columnindex 13
#define Db _p[14]
#define Db_columnindex 14
#define Db2 _p[15]
#define Db2_columnindex 15
#define v _p[16]
#define v_columnindex 16
#define _g _p[17]
#define _g_columnindex 17
#define _ion_ek	*_ppvar[0]._pval
#define _ion_ik	*_ppvar[1]._pval
#define _ion_dikdv	*_ppvar[2]._pval
 
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
 /* declaration of user functions */
 static void _hoc_a_tau(void);
 static void _hoc_a_inf(void);
 static void _hoc_b_tau(void);
 static void _hoc_b_inf(void);
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
 "setdata_kaDasoma", _hoc_setdata,
 "a_tau_kaDasoma", _hoc_a_tau,
 "a_inf_kaDasoma", _hoc_a_inf,
 "b_tau_kaDasoma", _hoc_b_tau,
 "b_inf_kaDasoma", _hoc_b_inf,
 "rates_kaDasoma", _hoc_rates,
 0, 0
};
#define a_tau a_tau_kaDasoma
#define a_inf a_inf_kaDasoma
#define b_tau b_tau_kaDasoma
#define b_inf b_inf_kaDasoma
 extern double a_tau( _threadargsprotocomma_ double );
 extern double a_inf( _threadargsprotocomma_ double );
 extern double b_tau( _threadargsprotocomma_ double );
 extern double b_inf( _threadargsprotocomma_ double );
 /* declare global and static user variables */
#define Vshift Vshift_kaDasoma
 double Vshift = -90;
#define Vmid_ina Vmid_ina_kaDasoma
 double Vmid_ina = -75;
#define Vmid_ac Vmid_ac_kaDasoma
 double Vmid_ac = -30;
#define h h_kaDasoma
 double h = 1;
#define k_ina k_ina_kaDasoma
 double k_ina = -7;
#define k_ac k_ac_kaDasoma
 double k_ac = 7;
#define m m_kaDasoma
 double m = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "Vmid_ac_kaDasoma", "mV",
 "k_ac_kaDasoma", "mV",
 "Vmid_ina_kaDasoma", "mV",
 "k_ina_kaDasoma", "mV",
 "Vshift_kaDasoma", "mV",
 "gbar_kaDasoma", "pS/microm2",
 "taurecov_kaDasoma", "ms",
 "i_kaDasoma", "mA/cm2",
 "g_kaDasoma", "pS/microm2",
 "atau_kaDasoma", "ms",
 "btau_kaDasoma", "ms",
 "ainf_kaDasoma", "1",
 "binf_kaDasoma", "1",
 0,0
};
 static double a0 = 0;
 static double b20 = 0;
 static double b0 = 0;
 static double delta_t = 0.01;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "Vmid_ac_kaDasoma", &Vmid_ac_kaDasoma,
 "k_ac_kaDasoma", &k_ac_kaDasoma,
 "Vmid_ina_kaDasoma", &Vmid_ina_kaDasoma,
 "k_ina_kaDasoma", &k_ina_kaDasoma,
 "Vshift_kaDasoma", &Vshift_kaDasoma,
 "m_kaDasoma", &m_kaDasoma,
 "h_kaDasoma", &h_kaDasoma,
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
 
#define _cvode_ieq _ppvar[3]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"kaDasoma",
 "gbar_kaDasoma",
 "taurecov_kaDasoma",
 0,
 "i_kaDasoma",
 "g_kaDasoma",
 "atau_kaDasoma",
 "btau_kaDasoma",
 "ainf_kaDasoma",
 "binf_kaDasoma",
 0,
 "a_kaDasoma",
 "b_kaDasoma",
 "b2_kaDasoma",
 0,
 0};
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 18, _prop);
 	/*initialize range parameters*/
 	gbar = 150;
 	taurecov = 25;
 	_prop->param = _p;
 	_prop->param_size = 18;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _KV4_DAsoma_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("k", -10000.);
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 18, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 kaDasoma /Users/cyprien/Desktop/EBrains_2025/Code/Mechanism/KV4_DAsoma.mod\n");
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
 static int _slist1[2], _dlist1[2];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {int _reset = 0; {
   rates ( _threadargscomma_ v ) ;
   Da = ( ainf - a ) / atau ;
   Db = ( binf - b ) / btau ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
 rates ( _threadargscomma_ v ) ;
 Da = Da  / (1. - dt*( ( ( ( - 1.0 ) ) ) / atau )) ;
 Db = Db  / (1. - dt*( ( ( ( - 1.0 ) ) ) / btau )) ;
  return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) { {
   rates ( _threadargscomma_ v ) ;
    a = a + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / atau)))*(- ( ( ( ainf ) ) / atau ) / ( ( ( ( - 1.0 ) ) ) / atau ) - a) ;
    b = b + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / btau)))*(- ( ( ( binf ) ) / btau ) / ( ( ( ( - 1.0 ) ) ) / btau ) - b) ;
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
 
double b_inf ( _threadargsprotocomma_ double _lV ) {
   double _lb_inf;
 _lb_inf = 1.0 / ( 1.0 + exp ( - ( _lV - Vmid_ina ) / k_ina ) ) ;
   
return _lb_inf;
 }
 
static void _hoc_b_inf(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  b_inf ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double a_tau ( _threadargsprotocomma_ double _lV ) {
   double _la_tau;
  _la_tau = 1.029 + ( 4.83 / ( 1.0 + exp ( ( _lV + 57.0 ) / 6.22 ) ) ) ;
    
return _la_tau;
 }
 
static void _hoc_a_tau(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  a_tau ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double b_tau ( _threadargsprotocomma_ double _lV ) {
   double _lb_tau;
  _lb_tau = taurecov + ( 120.0 + ( 78.4 / ( 1.0 + exp ( _lV + 68.5 ) / 5.95 ) ) - taurecov ) / ( 1.0 + exp ( ( - _lV + Vshift ) * 5.0 ) ) ;
    
return _lb_tau;
 }
 
static void _hoc_b_tau(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  b_tau ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
static int  rates ( _threadargsprotocomma_ double _lV ) {
   atau = a_tau ( _threadargscomma_ _lV ) ;
   ainf = a_inf ( _threadargscomma_ _lV ) ;
   btau = b_tau ( _threadargscomma_ _lV ) ;
   binf = b_inf ( _threadargscomma_ _lV ) ;
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
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ek = _ion_ek;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
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
  ek = _ion_ek;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_k_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 2, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
  a = a0;
  b2 = b20;
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
  ek = _ion_ek;
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   g = gbar * ( pow( a , m ) ) * ( pow( b , h ) ) ;
   i = ( 0.0001 ) * g * ( v - ek ) ;
   ik = i ;
   }
 _current += ik;

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
  ek = _ion_ek;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dik;
  _dik = ik;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ik += ik ;
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
  ek = _ion_ek;
 {   states(_p, _ppvar, _thread, _nt);
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = a_columnindex;  _dlist1[0] = Da_columnindex;
 _slist1[1] = b_columnindex;  _dlist1[1] = Db_columnindex;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/cyprien/Desktop/EBrains_2025/Code/Mechanism/KV4_DAsoma.mod";
static const char* nmodl_file_text = 
  "COMMENT \n"
  "\n"
  "Model for an kV4 cuRrent recorded in DA neurons.\n"
  "This current has two inactivation rates and a rapid rate for recovery from inactivation\n"
  "Activation and inactivation parameters are from amendola\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "\n"
  "NEURON\n"
  " {\n"
  "	 THREADSAFE\n"
  "  SUFFIX kaDasoma  USEION k READ ek WRITE ik \n"
  " \n"
  "RANGE gbar, g, i  \n"
  "\n"
  "RANGE atau, btau\n"
  "\n"
  "RANGE ainf, binf\n"
  "\n"
  "RANGE taurecov\n"
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
  "  gbar = 150 (pS/microm2)\n"
  "  Vmid_ac =-30 (mV)\n"
  "  k_ac = 7 (mV) \n"
  "  Vmid_ina = -75 (mV) \n"
  "  k_ina = -7 (mV)\n"
  "   taurecov=25 (ms)\n"
  "  Vshift=-90 (mV) : potential at which inactivation rates are replaced by taurecov\n"
  "   m=1         \n"
  "  h=1        : gate parameters according to the HH formalism (m*m*m*h)\n"
  "\n"
  "}\n"
  "\n"
  " ASSIGNED {\n"
  "  v	(mV)\n"
  "  ik 	(mA/cm2)\n"
  "  i 	(mA/cm2)\n"
  "  g	(pS/microm2)\n"
  "  atau (ms)\n"
  "  btau (ms)\n"
  "   \n"
  "  ainf (1)\n"
  "  binf (1)\n"
  " \n"
  "ek (mV)\n"
  " }\n"
  "\n"
  "\n"
  "STATE {a b b2}\n"
  "\n"
  "BREAKPOINT {\n"
  "  SOLVE states METHOD cnexp\n"
  "\n"
  " \n"
  "  g = gbar*(a^m)*(b^h)\n"
  "  i = (0.0001)*g*(v-ek)\n"
  "  ik = i\n"
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
  "  a' = (ainf-a)/atau\n"
  "  b' = (binf-b)/btau\n"
  "  \n"
  "}\n"
  "\n"
  "\n"
  "FUNCTION a_inf (V (mV)) () {\n"
  "\n"
  "  a_inf = 1/(1+exp(-(V-Vmid_ac)/k_ac))  : activation system (a*a*a)\n"
  "}\n"
  "\n"
  "FUNCTION b_inf (V (mV)) () {\n"
  "  b_inf = 1/(1+exp(-(V-Vmid_ina)/k_ina)) : inactivation system (b)\n"
  "}\n"
  "\n"
  "\n"
  "FUNCTION a_tau (V (mV)) (ms) {\n"
  "UNITSOFF\n"
  "\n"
  "a_tau= 1.029 + (4.83/(1+exp((V+57)/6.22)))\n"
  ": time constant of activation depends on V \n"
  "\n"
  "UNITSON\n"
  "}\n"
  "\n"
  "FUNCTION b_tau (V (mV)) (ms) {\n"
  "\n"
  "UNITSOFF\n"
  "\n"
  "  b_tau = taurecov + (120+(78.4/(1+exp(V+68.5)/5.95))-taurecov)/(1+exp((-V+Vshift)*5))\n"
  " : fast inactivation  \n"
  "\n"
  "\n"
  "UNITSON \n"
  "\n"
  "}\n"
  "\n"
  "\n"
  "\n"
  "PROCEDURE rates(V (mV)) {\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "atau=a_tau(V)\n"
  "\n"
  "ainf=a_inf(V)\n"
  "\n"
  "btau=b_tau(V)\n"
  "\n"
  "binf=b_inf(V)\n"
  "\n"
  "\n"
  "\n"
  "}\n"
  ;
#endif
