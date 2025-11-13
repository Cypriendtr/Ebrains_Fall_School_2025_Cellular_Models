TITLE calcium lVA channels for DA

COMMENT
 LOWthreshold calcium channel (L-type)

 
ENDCOMMENT

UNITS {
	(mM) = (milli/liter)
	(mV) = (millivolt)
	(mA) = (milliamp)
	(pS) =(picosiemens)
     FARADAY = (faraday) (coulomb)
           R = (k-mole) (joule/degC)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	THREADSAFE
	SUFFIX CAV13
	USEION ca READ eca WRITE ica
	RANGE  gbar,  iLCa
	
}

PARAMETER {
        v (mV)
	dt (ms)
       gbar  = 1.25 (pS/microm2)
			celsius
				  Vmid_ac =-31 (mV)
  k_ac = 7 (mV)
			iLCa=0
	
	
}

STATE {
        q
}

ASSIGNED { 
        ica (mA/cm2)
	qinf
		qtau (ms)
			eca (mV)
	
	
	
}

BREAKPOINT {
	LOCAL vghk
	SOLVE states METHOD cnexp
	:vghk = ghkg(v,cai,cao,2)
	iLCa = (gbar)*q*(v-eca)*0.0001
	ica  = iLCa
}


INITIAL {
	
        settables(v)
	q = qinf
	
}

DERIVATIVE states {  
	settables(v)  
	q' = (qinf-q)/qtau
	
}

PROCEDURE settables(v) {  :Computes rate and other constants at current v.
                          :Call once from HOC to initialize inf at resting v.
			  :Voltage shifts (for temp effects) of -8.25 and -14.67 added respt.
        TABLE qinf, qtau DEPEND celsius FROM -100 TO 100 WITH 400

                :"q" L Ca activation system
        qinf   = 1.0/(1.0 + exp((Vmid_ac - v)/ k_ac))
        qtau   =1/((-0.20876*(v+39.26)/(exp(-(v+39.26)/4.111)-1)+(0.9444*exp(-(v+15.38)/224.1))))
					
      
}



: INCLUDE "ghk.inc"


