const e20 = 2.2; // rate of change of saturated VP with T at 20C
const maxSoilCond = 0.00250;
const lambda = 2460000.0; // latent heat of vapourisation of H2O (J/kg)
const rhoAir = 1.2; // density of air, kg/m3
const VPDconv = 0.000622; // convert VPD to saturation deficit = 18/29/1000

const dayOfYear = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349];
const daysInMonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
