#pragma once
#include "HepMC3/GenEvent.h"
#include "HepMC3/Attribute.h"
#include "HepMC3/GenParticle.h"
#include <HepMC3/GenVertex.h>
#include "EventInitilizers.h"


class HepMCReader {
public:
  HepMCReader() = default;
  ~HepMCReader();

  // Reads the next event and fills evtIn. Returns false on EOF or failure.
  void ReadAndFillEvent(EventInitilizers& evtIn, HepMC3::GenEvent& evt);

private:
  template <class T>
  inline void assign4(std::array<double,4>& dst, const std::vector<T>& src) {
    if (src.size() < 4) throw std::runtime_error("attribute has <4 elements");
    dst[0] = static_cast<double>(src[3]);
    dst[1] = static_cast<double>(src[0]);
    dst[2] = static_cast<double>(src[1]);
    dst[3] = static_cast<double>(src[2]);
  }

  // helper: copy HepMC3 4-vector into your {E,px,py,pz} array
  static inline void set4(std::array<double,4>& dst, const HepMC3::FourVector& v) {
    dst[0] = v.e();  dst[1] = v.px();  dst[2] = v.py();  dst[3] = v.pz();
  }

};