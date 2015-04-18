#ifndef sph_writer_hgztek43_hh
#define sph_writer_hgztek43_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

template<class PartikelArray>
struct Writer {

  typedef typename PartikelArray::value_type particle_type;

  enum { 
    SAVE_MOVING = 1, 
    SAVE_BOUNDARY = 2,
    SAVE_OUT = 4
  };

  void findParticlesToSave(PartikelArray const &particle, 
                           long const saveFlags,
                           std::valarray<unsigned> &pindex2) const {
    size_t const nparticle = particle.size();
    size_t pToSaveCount = 0;
    {
      std::valarray<unsigned> pindex(nparticle);
      for(size_t i = 0; i < nparticle; ++i) {
        particle_type const &p = particle[i];
        if ((p.isBoundary() and (saveFlags & SAVE_BOUNDARY))
            or (!p.isBoundary() and (saveFlags & SAVE_MOVING))) {
          if ((p.isOut() and (saveFlags & SAVE_OUT))
              or not p.isOut()) {
            pindex[pToSaveCount] = i;
            ++pToSaveCount;
          }
        }
      }
      pindex2.resize(pToSaveCount);
      pindex2 = pindex[std::slice(0, pToSaveCount, 1)];
    }
  }

};

#endif
