#ifndef ANALYSIS_UTILS_H
#define ANALYSIS_UTILS_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <map>

#define PI 3.14159265359



float DeltaEta(const float& eta1, const float& eta2);
float DeltaPhi(const float& phi1, const float& phi2);
float DeltaR(const float& eta1, const float& phi1,
             const float& eta2, const float& phi2);

#endif
