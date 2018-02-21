#ifndef SRC_ISS_HADRON_H_
#define SRC_ISS_HADRON_H_

struct iSS_Hadron {
     int pid;
     double mass;
     double E, px, py, pz;
     double t, x, y, z;
};

#endif  // SRC_ISS_HADRON_H_
