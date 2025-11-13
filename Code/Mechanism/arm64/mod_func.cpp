#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
#if defined(__cplusplus)
extern "C" {
#endif

extern void _cad_reg(void);
extern void _CAV13_reg(void);
extern void _IHnew_reg(void);
extern void _kca_callaway_reg(void);
extern void _KV4_DA_reg(void);
extern void _KV4_DAsoma_reg(void);
extern void _KVDRNTS_reg(void);
extern void _NA12_reg(void);
extern void _passive_reg(void);

void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"cad.mod\"");
    fprintf(stderr, " \"CAV13.mod\"");
    fprintf(stderr, " \"IHnew.mod\"");
    fprintf(stderr, " \"kca_callaway.mod\"");
    fprintf(stderr, " \"KV4_DA.mod\"");
    fprintf(stderr, " \"KV4_DAsoma.mod\"");
    fprintf(stderr, " \"KVDRNTS.mod\"");
    fprintf(stderr, " \"NA12.mod\"");
    fprintf(stderr, " \"passive.mod\"");
    fprintf(stderr, "\n");
  }
  _cad_reg();
  _CAV13_reg();
  _IHnew_reg();
  _kca_callaway_reg();
  _KV4_DA_reg();
  _KV4_DAsoma_reg();
  _KVDRNTS_reg();
  _NA12_reg();
  _passive_reg();
}

#if defined(__cplusplus)
}
#endif
