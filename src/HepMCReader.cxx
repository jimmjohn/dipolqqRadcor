#include "HepMCReader.h"

HepMCReader::~HepMCReader() {
    // Destructor implementation
}

void HepMCReader::ReadAndFillEvent(EventInitilizers& evtIn, HepMC3::GenEvent& evt) {
    // Implementation of reading from HepMC3 event and filling evtIn

  if (auto a = evt.attribute<HepMC3::FloatAttribute>("SpinWT"))
        evtIn.wt = a->value();

  if (auto a = evt.attribute<HepMC3::FloatAttribute>("SpinWThelApprox"))
        evtIn.wt_approx = a->value();

  // Reading approximate helicities
  for (auto& p : evt.particles()) {
     if (!p) continue;
        if (auto a = p->attribute<HepMC3::IntAttribute>("approximateHelicity")) {
            if(evtIn.heTaum==0) evtIn.heTaum = a->value();
            else evtIn.heTaul = a->value();
        }
  }

  //Reading Polarimetic vectors in lab frame
  for (auto& p : evt.particles()) {
      if (!p) continue;
      if (auto a = p->attribute<HepMC3::VectorFloatAttribute>("polarimetricInLabFrame")) {
          if(evtIn.pvTaum_found==false) {
              assign4(evtIn.m_HvCloneTaum, a->value());
              evtIn.pvTaum_found = true;
          }
          else {
              assign4(evtIn.m_HvCloneTaul, a->value());
              evtIn.pvTaul_found = true;
          }
      }
  }

  //Reading the tau particle and its daughters

  for(auto& p : evt.particles()) {
      if(p->pdg_id() == 15) {
          set4(evtIn.P1, p->momentum());
          auto vtx = p->end_vertex();
          for (auto& d : vtx->particles_out()) {
              if (!d) continue;
              if(d->pdg_id() == -211) set4(evtIn.PIPM, d->momentum());
              else if(d->pdg_id() == 111) set4(evtIn.PIZM, d->momentum());
          }
      }
      else if(p->pdg_id() == -15) {
          set4(evtIn.P2, p->momentum());
          auto vtx = p->end_vertex();
          for (auto& d : vtx->particles_out()) {
              if (!d) continue;
              if(d->pdg_id() == 211) set4(evtIn.PIPL, d->momentum());
              else if(d->pdg_id() == 111) set4(evtIn.PIZL, d->momentum());
          }
      }
      else if(p->pdg_id() == 11) {
            evtIn.beamEnergy += p->momentum().e();
      } else if(p->pdg_id() == -11) {
            evtIn.beamEnergy += p->momentum().e();
      }
  }
}
