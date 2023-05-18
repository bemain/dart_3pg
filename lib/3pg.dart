import 'dart:math';

import 'constrants.dart';

class SiteData {
  SiteData(
    this.lat,
    this.altitude,
    this.soilClass,
    this.asw,
    this.asw_min,
    this.asw_max,
    this.year_i,
    this.month_i,
  );

  /// site latitude
  final double lat;

  /// altitude of the site location, m
  final int altitude;

  /// soil parameters for soil class
  final int soilClass;

  /// available soil water
  final double asw;

  /// minimum available soil water
  final double asw_min;

  /// maximum available soil water
  final double asw_max;

  /// initial year when the simulations starts
  final int year_i;

  /// initial month when the simulation starts
  final int month_i;
}

class SpeciesData {
  SpeciesData(
    this.year_p,
    this.month_p,
    this.fertility,
    this.stems_n_i,
    this.biom_stem_i,
    this.biom_root_i,
    this.biom_foliage_i,
  );

  /// year when species was planted
  final List<int> year_p;

  /// month when species was planted
  final List<int> month_p;

  /// initial site fertility rating for a given species
  final List<double> fertility;

  /// initial stand stocking for a species
  final List<double> stems_n_i;

  /// initial stem biomass for a species
  final List<double> biom_stem_i;

  /// initial root biomass for a species
  final List<double> biom_root_i;

  /// initial foliage biomass for a species
  final List<double> biom_foliage_i;
}

class ClimateData {
  ClimateData(
    this.tmp_min,
    this.tmp_max,
    this.tmp_ave,
    this.prcp,
    this.solar_rad,
    this.frost_days,
    this.vpd_day,
    this.co2,
    this.d13catm,
  );

  /// minimum daily temperature
  final double tmp_min;

  /// maximum daily temperature
  final double tmp_max;

  /// average daily temperature
  final double tmp_ave;

  /// monthly precipitation sum
  final double prcp;

  /// mean daily incident solar radiation
  final double solar_rad;

  /// number of frost days per month
  final double frost_days;

  final double vpd_day;

  /// atmospheric CO2
  final double co2;

  /// added d13C of atmospheric CO2 (per mil)
  final double d13catm;
}

class Settings {
  Settings(
    this.light_model,
    this.transp_model,
    this.phys_model,
    this.height_model,
    this.correct_bias,
    this.calculate_d13c,
  );

  /// 1 - 3PGpjs; 2 - 3PGmix
  final int light_model;

  /// 1 - 3PGpjs; 2 - 3PGmix
  final int transp_model;

  /// 1 - 3PGpjs; 2 - 3PGmix
  final int phys_model;

  /// 1 - linear; 2 - non-linear
  final int height_model;

  /// 0 - no; 1 - 3PGmix
  final int correct_bias;

  /// 0 - no; 1 - 3PGmix
  final int calculate_d13c;
}

class ParsI {
  /// One value for each species.
  ParsI(
      this.pFS2,
      this.pFS20,
      this.aWs,
      this.nWs,
      this.pRx,
      this.pRn,
      this.gammaF1,
      this.gammaF0,
      this.tgammaF,
      this.gammaR,
      this.leafgrow,
      this.leaffall,
      this.Tmin,
      this.Topt,
      this.Tmax,
      this.kF,
      this.SWconst0,
      this.SWpower0,
      this.fCalpha700,
      this.fCg700,
      this.m0,
      this.fN0,
      this.fNn,
      this.MaxAge,
      this.nAge,
      this.rAge,
      this.gammaN1,
      this.gammaN0,
      this.tgammaN,
      this.ngammaN,
      this.wSx1000,
      this.thinPower,
      this.mF,
      this.mR,
      this.mS,
      this.SLA0,
      this.SLA1,
      this.tSLA,
      this.k,
      this.fullCanAge,
      this.MaxIntcptn,
      this.LAImaxIntcptn,
      this.cVPD,
      this.alphaCx,
      this.y,
      this.MinCond,
      this.MaxCond,
      this.LAIgcx,
      this.CoeffCond,
      this.BLcond,
      this.RGcGW,
      this.D13CTissueDif,
      this.aFracDiffu,
      this.bFracRubi,
      this.fracBB0,
      this.fracBB1,
      this.tBB,
      this.rho0,
      this.rho1,
      this.tRho,
      this.CrownShape,
      this.aH,
      this.nHB,
      this.nHC,
      this.aV,
      this.nVB,
      this.nVH,
      this.nVBH,
      this.aK,
      this.nKB,
      this.nKH,
      this.nKC,
      this.nKrh,
      this.aHL,
      this.nHLB,
      this.nHLL,
      this.nHLC,
      this.nHLrh,
      this.Qa,
      this.Qb,
      this.gDM_mol,
      this.molPAR_MJ);

  final List<double> pFS2;
  final List<double> pFS20;
  final List<double> aWs;
  final List<double> nWs;
  final List<double> pRx;
  final List<double> pRn;
  final List<double> gammaF1;
  final List<double> gammaF0;
  final List<double> tgammaF;
  final List<double> gammaR;
  final List<int> leafgrow;
  final List<int> leaffall;
  final List<double> Tmin;
  final List<double> Topt;
  final List<double> Tmax;
  final List<double> kF;
  final List<double> SWconst0;
  final List<double> SWpower0;
  final List<double> fCalpha700;
  final List<double> fCg700;
  final List<double> m0;
  final List<double> fN0;
  final List<double> fNn;
  final List<double> MaxAge;
  final List<double> nAge;
  final List<double> rAge;
  final List<double> gammaN1;
  final List<double> gammaN0;
  final List<double> tgammaN;
  final List<double> ngammaN;
  final List<double> wSx1000;
  final List<double> thinPower;
  final List<double> mF;
  final List<double> mR;
  final List<double> mS;
  final List<double> SLA0;
  final List<double> SLA1;
  final List<double> tSLA;
  final List<double> k;
  final List<double> fullCanAge;
  final List<double> MaxIntcptn;
  final List<double> LAImaxIntcptn;
  final List<double> cVPD;
  final List<double> alphaCx;
  final List<double> y;
  final List<double> MinCond;
  final List<double> MaxCond;
  final List<double> LAIgcx;
  final List<double> CoeffCond;
  final List<double> BLcond;
  final List<double> RGcGW;
  final List<double> D13CTissueDif;
  final List<double> aFracDiffu;
  final List<double> bFracRubi;
  final List<double> fracBB0;
  final List<double> fracBB1;
  final List<double> tBB;
  final List<double> rho0;
  final List<double> rho1;
  final List<double> tRho;
  final List<int> CrownShape;
  final List<double> aH;
  final List<double> nHB;
  final List<double> nHC;
  final List<double> aV;
  final List<double> nVB;
  final List<double> nVH;
  final List<double> nVBH;
  final List<double> aK;
  final List<double> nKB;
  final List<double> nKH;
  final List<double> nKC;
  final List<double> nKrh;
  final List<double> aHL;
  final List<double> nHLB;
  final List<double> nHLL;
  final List<double> nHLC;
  final List<double> nHLrh;
  final List<double> Qa;
  final List<double> Qb;
  final List<double> gDM_mol;
  final List<double> molPAR_MJ;
}

class ParsB {
  ParsB(
      this.Dscale0,
      this.DscaleB,
      this.Dscalerh,
      this.Dscalet,
      this.DscaleC,
      this.Dshape0,
      this.DshapeB,
      this.Dshaperh,
      this.Dshapet,
      this.DshapeC,
      this.Dlocation0,
      this.DlocationB,
      this.Dlocationrh,
      this.Dlocationt,
      this.DlocationC,
      this.wsscale0,
      this.wsscaleB,
      this.wsscalerh,
      this.wsscalet,
      this.wsscaleC,
      this.wsshape0,
      this.wsshapeB,
      this.wsshaperh,
      this.wsshapet,
      this.wsshapeC,
      this.wslocation0,
      this.wslocationB,
      this.wslocationrh,
      this.wslocationt,
      this.wslocationC);

  final List<double> Dscale0;
  final List<double> DscaleB;
  final List<double> Dscalerh;
  final List<double> Dscalet;
  final List<double> DscaleC;

  final List<double> Dshape0;
  final List<double> DshapeB;
  final List<double> Dshaperh;
  final List<double> Dshapet;
  final List<double> DshapeC;

  final List<double> Dlocation0;
  final List<double> DlocationB;
  final List<double> Dlocationrh;
  final List<double> Dlocationt;
  final List<double> DlocationC;

  final List<double> wsscale0;
  final List<double> wsscaleB;
  final List<double> wsscalerh;
  final List<double> wsscalet;
  final List<double> wsscaleC;

  final List<double> wsshape0;
  final List<double> wsshapeB;
  final List<double> wsshaperh;
  final List<double> wsshapet;
  final List<double> wsshapeC;

  final List<double> wslocation0;
  final List<double> wslocationB;
  final List<double> wslocationrh;
  final List<double> wslocationt;
  final List<double> wslocationC;
}

class ParsS {
  ParsS(
    this.aH,
    this.nHB,
    this.nHC,
    this.aV,
    this.nVB,
    this.nVH,
    this.nVBH,
    this.aK,
    this.nKB,
    this.nKH,
    this.nKC,
    this.nKrh,
    this.aHL,
    this.nHLB,
    this.nHLL,
    this.nHLC,
    this.nHLrh,
  );

  final List<double> aH;
  final List<double> nHB;
  final List<double> nHC;
  final List<double> aV;
  final List<double> nVB;
  final List<double> nVH;
  final List<double> nVBH;
  final List<double> aK;
  final List<double> nKB;
  final List<double> nKH;
  final List<double> nKC;
  final List<double> nKrh;
  final List<double> aHL;
  final List<double> nHLB;
  final List<double> nHLL;
  final List<double> nHLC;
  final List<double> nHLrh;
}

/// Age of each species and month
List<List<double>> age = [];

/// Age of each species used for calculating modifiers (one month less than s_age)
List<List<double>> age_m = [];

List<double> stems_n = [];

/// potential number of stems per ha
List<double> stems_n_ha = [];

/// stand level basal area
List<double> basal_area = [];

/// proportion of basal area
List<double> basal_area_prop = [];

/// average tree DBH, cm
List<double> dbh = [];

/// average tree height, m
List<double> height = [];

double Height_max = 0;

/// ***DF mean live-crown length (m) of a species
List<double> crown_length = [];

/// ***DF mean crown diameter (m)
List<double> crown_width = [];

List<double> volume = [];
List<double> volume_mai = [];
List<double> volume_old = [];
List<double> volume_cum = [];
List<double> volume_change = [];

List<double> competition_total = [];

/// Specific leaf area
List<List<double>> SLA = [];

/// Fraction of stem biomass as branch and bark
List<List<double>> fracBB = [];

/// Whole-tree basic density
List<List<double>> wood_density = [];

/// Canopy LAI (mean annual LAI if output time step is annual, and final year LAI if step is whole rotation)
List<double> LAI = [];

/// total competition of the forest
List<double> lai_total = [];

/// species specific proportion of lai
List<double> LAI_per = [];

/// leaf area above the given species
List<double> lai_above = [];

List<double> canopy_vol_frac = [];

/// the ratio of mean tree leaf area (m2) to crownSA (m2)
List<double> lai_sa_ratio = [];

List<int> layer_id = [];

/// Foliage biomass
List<double> biom_foliage = [];

List<double> biom_foliage_debt = [];

/// Root biomass
List<double> biom_root = [];

/// Stem biomass, including branches and bark
List<double> biom_stem = [];

/// average tree stem mass
List<double> biom_tree = [];

/// Max. mean tree stem mass at current stocking
List<double> biom_tree_max = [];

List<double> biom_incr_foliage = [];
List<double> biom_incr_root = [];
List<double> biom_incr_stem = [];

/// Litter fall
List<double> biom_loss_foliage = [];
List<double> biom_loss_root = [];

dynamic s_3PG_f(
  SiteData siteInputs,
  List<SpeciesData> speciesInputs,
  List<ClimateData> forcingInputs,
  List<List<List<double>>> managementInputs,
  ParsI pars_i,
  ParsB pars_b,
  List<int> t_t,
  Settings settings, // settings for the models
) {
  int n_m = forcingInputs.length; // number of months
  int n_sp = speciesInputs.length; // number of species
  int n_man = managementInputs.length; // number of management interactions

  // ====== INITIALIZATION ======
}

bool f_dormant(int month, int leafgrow, int leaffall) {
  // This is called if the leafgrow parameter is not 0, and hence the species is Deciduous
  // This is true if "currentmonth" is part of the dormant season
  if (leafgrow > leaffall) {
    // check which hemisphere
    if (month >= leaffall && month <= leafgrow) {
      // growing at winter
      return true;
    } else if (leafgrow < leaffall) {
      if (month < leafgrow || month >= leaffall) {
        // growing at summer
        return true;
      }
    }
  }

  return false;
}

List<double> f_exp(
  int n_m, // Not used in the method
  List<double> x,
  double g0,
  double gx,
  double tg,
  double ng,
) {
  if (tg != 0) {
    return x
        .map((val) => gx + (g0 - gx) * exp(-ln2 * pow(val / tg, ng)))
        .toList();
  }

  return List.filled(x.length, gx);
}

List<double> f_exp_foliage(
  int n_m, // Not used in the method
  List<double> x,
  double f1,
  double f0,
  double tg,
) {
  if (tg * f1 == 0) {
    return List.filled(x.length, f1);
  } else {
    double kg = 12.0 * log(1.0 + f1 / f0) / tg;
    return x.map((val) => f1 * f0 / (f0 + (f1 - f0) * exp(-kg * val))).toList();
  }
}

List<double> f_gamma_dist(List<double> x, int n) {
  return x
      .map((val) =>
          pow(val, (val - 0.5)) *
          exp(-val) *
          sqrt(2.0 * pi) *
          (1.0 +
                  1.0 / (12.0 * val) +
                  1.0 / (288.0 * pow(val, 2.0)) -
                  139.0 / (51840.0 * pow(val, 3.0)) -
                  571.0 / (2488320.0 * pow(val, 4.0)))
              .toDouble())
      .toList();
}

List<double> f_get_daylength(double Lat) {
  List<double> day_length = List<double>.filled(12, 0);
  List<double> sinDec = List<double>.filled(12, 0);
  List<double> cosH0 = List<double>.filled(12, 0);

  double SLAt = sin(pi * Lat / 180);
  double cLat = cos(pi * Lat / 180);

  for (int i = 0; i < 12; i++) {
    sinDec[i] = 0.4 * sin(0.0172 * (dayOfYear[i] - 80));
    cosH0[i] = -sinDec[i] * SLAt / (cLat * sqrt(1 - pow(sinDec[i], 2)));
    day_length[i] = acos(cosH0[i]) / pi;

    if (cosH0[i] > 1) day_length[i] = 0;
    if (cosH0[i] < -1) day_length[i] = 1;
  }

  return day_length;
}

/// function to allocate each tree to the layer based on height and crown heigh
/// First layer (1) is the highest
/// According to Forrester, D.I., Guisasola, R., Tang, X. et al. For. Ecosyst. (2014) 1: 17.
/// Calculations based on example https://it.mathworks.com/matlabcentral/answers/366626-overlapping-time-intervals
List<int> f_get_layer(int n_sp, List<double> height, List<double> Heightcrown) {
  List<int> layer_id = List<int>.filled(n_sp, 0);
  List<double> Height_all = List<double>.filled(n_sp * 2, 0.0);
  List<int> Height_ind = List<int>.filled(n_sp * 2, 0);
  List<int> ones = List<int>.filled(n_sp * 2, -1);
  List<int> ones_sum = List<int>.filled(n_sp * 2, 0);
  List<double> Height_layer;
  int i;
  int n_l;

  // Sort all height and crown height
  for (i = 0; i < n_sp; i++) {
    Height_all[i] = Heightcrown[i];
    Height_all[i + n_sp] = height[i];
  }

  Height_ind = f_orderId(Height_all);

  // Assign index order for further calculations
  for (i = 0; i < n_sp * 2; i++) {
    ones[i] = 1;
  }
  ones = ones.sublist(0, Height_ind.length);

  // Cumulative sum
  ones_sum[0] = ones[0];
  for (i = 1; i < n_sp * 2; i++) {
    ones_sum[i] = ones_sum[i - 1] + ones[i];
  }

  // Max height of each layer
  n_l = ones_sum.where((element) => element == 0).length;
  Height_layer = List<double>.filled(n_l, 0.0);
  for (i = 0; i < n_l; i++) {
    Height_layer[i] = Height_all[Height_ind[i]];
  }

  // Assign layer to each species
  for (i = 0; i < n_sp; i++) {
    layer_id[i] = 1;
  }
  if (n_l > 1) {
    for (i = 0; i < n_l - 1; i++) {
      for (int j = 0; j < n_sp; j++) {
        if (height[j] > Height_layer[i]) {
          layer_id[j] = i + 2;
        }
      }
    }
  }

  // Revert the order, so highest trees are 1 layer and lowest is n
  int max_layer = layer_id.reduce(max);
  for (i = 0; i < n_sp; i++) {
    layer_id[i] = max_layer - layer_id[i] + 1;
  }

  return layer_id;
}

List<double> f_get_layer_sum(
    int n_sp, int nLayers, List<double> x, List<int> layer_id) {
  List<double> y = List<double>.filled(n_sp, 0.0);
  int i;

  for (i = 0; i < nLayers; i++) {
    for (int j = 0; j < n_sp; j++) {
      if (layer_id[j] == i + 1) {
        y[j] += x[j];
      }
    }
  }

  return y;
}

double f_get_mortality(
    double stems_n, double WS, double mS, double wSx1000, double thinPower) {
  double mort_n;
  double accuracy = 1.0 / 1000.0;
  int i;
  double fN, dfN, dN, n, x1, x2;

  n = stems_n / 1000.0;
  x1 = 1000.0 * mS * WS / stems_n;
  i = 0;

  do {
    i = i + 1;

    if (n <= 0.0) {
      break; // added in 3PG+
    }

    x2 = wSx1000 * pow(n, (1.0 - thinPower));
    fN = x2 - x1 * n - (1.0 - mS) * WS;
    dfN = (1.0 - thinPower) * x2 / n - x1;
    dN = -fN / dfN;
    n = n + dN;

    if (dN.abs() <= accuracy || i >= 5) {
      break;
    }
  } while (true);

  mort_n = stems_n - 1000.0 * n;

  return mort_n;
}

List<double> f_get_solarangle(List<double> dayOfYear, double Lat) {
  List<double> solarangle = List<double>.filled(12, 0.0);
  double secondxaxisintercept, firstxaxisintercept;
  List<double> gamma = List<double>.filled(12, 0.0);
  List<double> declinationangle = List<double>.filled(12, 0.0);
  List<double> szaprep = List<double>.filled(12, 0.0);
  List<double> solarzenithangle = List<double>.filled(12, 0.0);

  secondxaxisintercept =
      0.0018 * pow(Lat, 3.0) - 0.0031 * pow(Lat, 2.0) + 2.3826 * Lat + 266.62;
  firstxaxisintercept =
      -0.0018 * pow(Lat, 3.0) + 0.0021 * pow(Lat, 2.0) - 2.3459 * Lat + 80.097;

  gamma =
      List<double>.generate(12, (i) => 2.0 * pi / 365.0 * (dayOfYear[i] - 1));

  for (int i = 0; i < 12; i++) {
    declinationangle[i] = 0.006918 -
        (0.399912 * cos(gamma[i])) +
        0.070257 * sin(gamma[i]) -
        0.006758 * cos(2.0 * gamma[i]) +
        0.000907 * sin(2.0 * gamma[i]) -
        0.002697 * cos(3.0 * gamma[i]) +
        0.00148 * sin(3.0 * gamma[i]);

    szaprep[i] = sin(pi / 180.0 * Lat * (-1.0)) * sin(declinationangle[i]) +
        cos(pi / 180.0 * Lat * (-1.0)) * cos(declinationangle[i]);
    solarzenithangle[i] = 180.0 /
        pi *
        (atan(-szaprep[i] / sqrt(-szaprep[i] * szaprep[i] + 1.0)) +
            2.0 * atan(1.0));
  }

  solarangle = solarzenithangle;

  if (Lat >= 0.0 && Lat <= 23.4) {
    for (int i = 0; i < 12; i++) {
      if (dayOfYear[i] > secondxaxisintercept ||
          dayOfYear[i] < firstxaxisintercept) {
        solarangle[i] = -1.0 * solarzenithangle[i];
      }
    }
  }

  if (Lat >= -23.4 && Lat < 0.0) {
    for (int i = 0; i < 12; i++) {
      if (dayOfYear[i] > firstxaxisintercept &&
          dayOfYear[i] < secondxaxisintercept) {
        solarangle[i] = -1.0 * solarzenithangle[i];
      }
    }
  }

  return solarangle;
}

List<int> f_orderId(List<double> x) {
  List<int> id = List<int>.filled(x.length, 0);
  int i, n, imin, temp1;
  double temp2;
  List<double> x2 = List<double>.from(x);

  n = x.length;

  for (i = 0; i < n; i++) {
    id[i] = i + 1;
  }

  for (i = 0; i < n - 1; i++) {
    imin =
        x2.sublist(i).indexOf(x2.sublist(i).reduce((a, b) => a < b ? a : b)) +
            i;
    if (imin != i) {
      temp2 = x2[i];
      x2[i] = x2[imin];
      x2[imin] = temp2;

      temp1 = id[i];
      id[i] = id[imin];
      id[imin] = temp1;
    }
  }

  return id;
}

List<double> p_min_max(List<double> x, double mn, double mx) {
  int n = x.length;
  List<double> out = List<double>.filled(n, 0);

  for (int i = 0; i < n; i++) {
    if (x[i] > mx) {
      out[i] = mx;
    } else if (x[i] < mn) {
      out[i] = mn;
    } else {
      out[i] = x[i];
    }
  }

  return out;
}

/// Returns canopy_cover and apar in a list
List<List<double>> s_light_3pgpjs(
  int n_sp,
  List<double> age,
  List<double> fullCanAge,
  List<double> k,
  List<double> lai,
  List<double> solar_rad,
  List<int> days_in_month,
) {
  List<double> canopy_cover = List<double>.filled(n_sp, 0.0);
  List<double> apar = List<double>.filled(n_sp, 0.0);

  for (int i = 0; i < n_sp; i++) {
    canopy_cover[i] = 1.0;
    if (fullCanAge[i] > 0.0 && age[i] < fullCanAge[i]) {
      canopy_cover[i] = (age[i] + 0.01) / fullCanAge[i];
    }
  }

  for (int i = 0; i < n_sp; i++) {
    double lightIntcptn = 1.0 - exp(-k[i] * lai[i] / canopy_cover[i]);
    double RADt = solar_rad[i] * days_in_month[i];
    apar[i] = RADt * lightIntcptn * canopy_cover[i];
  }

  return [canopy_cover, apar];
}

/// Returns apar, lai_above, fi, lambda_v, lambda_h, canopy_vol_frac, layer_id, lai_sa_ratio in a list
List s_light_3pgmix(
  int n_sp,
  List<double> height,
  List<double> crown_length,
  List<double> crown_width,
  List<double> lai,
  List<double> stems_n,
  double solar_rad,
  List<int> CrownShape,
  List<double> k,
  double solarAngle,
  int days_in_month,
) {
  List<double> Heightmidcrown = List<double>.filled(n_sp, 0.0);
  List<double> Heightcrown = List<double>.filled(n_sp, 0.0);
  List<double> CrownSA = List<double>.filled(n_sp, 0.0);
  List<double> Crownvolume = List<double>.filled(n_sp, 0.0);
  List<double> Height_max_l = List<double>.filled(n_sp, 0.0);
  List<double> Heightcrown_min_l = List<double>.filled(n_sp, 0.0);
  List<double> Heightmidcrown_l = List<double>.filled(n_sp, 0.0);
  List<double> Heightmidcrown_r = List<double>.filled(n_sp, 0.0);
  List<double> kL_l = List<double>.filled(n_sp, 0.0);
  List<double> lambdaV_l = List<double>.filled(n_sp, 0.0);
  List<double> kLSweightedave = List<double>.filled(n_sp, 0.0);
  List<double> aparl = List<double>.filled(n_sp, 0.0);
  double RADt = 0.0;
  List<double> LAI_l = List<double>.filled(n_sp, 0.0);

  List<double> apar = List<double>.filled(n_sp, 0.0);
  List<double> lai_above = List<double>.filled(n_sp, 0.0);
  List<double> fi = List<double>.filled(n_sp, 0.0);
  List<double> lambda_v = List<double>.filled(n_sp, 0.0);
  List<double> lambda_h = List<double>.filled(n_sp, 0.0);
  List<double> canopy_vol_frac = List<double>.filled(n_sp, 0.0);
  List<int> layer_id = List<int>.filled(n_sp, 0);
  List<double> lai_sa_ratio = List<double>.filled(n_sp, 0.0);

  for (int i = 0; i < n_sp; i++) {
    Heightcrown[i] = height[i] - crown_length[i];
    Heightmidcrown[i] = height[i] - crown_length[i] / 2;
  }

  /// Calculate the crown area and volume
  /// We only do it for species that have LAI, otherwise it stays 0 as was initialized above
  for (int i = 0; i < n_sp; i++) {
    if (lai[i] > 0.0) {
      if (CrownShape[i] == 1) {
        // Cone shaped
        CrownSA[i] = pi * pow((crown_width[i] / 2.0), 2) +
            pi *
                crown_width[i] /
                2.0 *
                pow(pow((crown_width[i] / 2.0), 2) + pow(crown_length[i], 2),
                    0.5);
        Crownvolume[i] =
            (pi * crown_width[i] * crown_width[i] * crown_length[i]) / 12.0;
      } else if (CrownShape[i] == 2) {
        // Ellipsoid
        CrownSA[i] = 4.0 *
            pi *
            pow(
                (((pow((crown_width[i] / 2.0), 1.6075) *
                            pow((crown_width[i] / 2.0), 1.6075)) +
                        (pow((crown_width[i] / 2.0), 1.6075) *
                            pow((crown_length[i] / 2.0), 1.6075)) +
                        (pow((crown_width[i] / 2.0), 1.6075) *
                            pow((crown_length[i] / 2.0), 1.6075))) /
                    3.0),
                (1.0 / 1.6075));
        Crownvolume[i] =
            (pi * crown_width[i] * crown_width[i] * crown_length[i] * 4.0) /
                24.0;
      } else if (CrownShape[i] == 3) {
        // Half-ellipsoid
        CrownSA[i] = pi * pow((crown_width[i] / 2.0), 2) +
            (4.0 *
                    pi *
                    pow(
                        (((pow((crown_width[i] / 2.0), 1.6075) *
                                    pow((crown_width[i] / 2.0), 1.6075)) +
                                (pow((crown_width[i] / 2.0), 1.6075) *
                                    pow(crown_length[i], 1.6075)) +
                                (pow((crown_width[i] / 2.0), 1.6075) *
                                    pow(crown_length[i], 1.6075))) /
                            3.0),
                        (1.0 / 1.6075))) /
                2.0;
        Crownvolume[i] =
            (pi * crown_width[i] * crown_width[i] * crown_length[i] * 4.0) /
                24.0;
      } else if (CrownShape[i] == 4) {
        // Rectangular
        CrownSA[i] = crown_width[i] * crown_width[i] * 2.0 +
            crown_width[i] * crown_length[i] * 4.0;
        Crownvolume[i] = crown_width[i] * crown_width[i] * crown_length[i];
      }
    }
  }

  // Calculate the ratio of tree leaf area to crown surface area restrict kLS to 1
  for (int i = 0; i < n_sp; i++) {
    lai_sa_ratio[i] = lai[i] * 10000.0 / stems_n[i] / CrownSA[i];
    if (lai[i] == 0.0) {
      lai_sa_ratio[i] = 0.0;
    }
  }

// Separate trees into layers
  layer_id = f_get_layer(n_sp, height, Heightcrown);
  int nLayers = layer_id.reduce((max, value) => max > value ? max : value);

//Now calculate the proportion of the canopy space that is filled by the crowns. The canopy space is the
  // volume between the top and bottom of a layer that is filled by crowns in that layer.
  // We calculate it only for the trees that have LAI and are in that particular year. Thus the tree can be in that
  // layer, but currently will not have LAI
  for (int i = 1; i <= nLayers; i++) {
    double max = 0;
    for (int j = 0; j < height.length; j++) {
      if (layer_id[j] == i && lai[j] != 0.0 && height[j] > max) max = height[j];
    }
    double Height_max_l = max;

    double min = 0;
    for (int j = 0; j < Heightcrown.length; j++) {
      if (layer_id[j] == i && lai[j] != 0.0 && Heightcrown[j] < max)
        max = Heightcrown[j];
    }
    double Heightcrown_min_l = min;
  }

// Calculate the ratio between the mid height of the given species and the mid height of the layer
  for (int i = 0; i < n_sp; i++) {
    // sum the canopy volume fraction per layer and save it at each species
    canopy_vol_frac[i] = Crownvolume[i] *
        stems_n[i] /
        ((Height_max_l[i] - Heightcrown_min_l[i]) * 10000);
    canopy_vol_frac = f_get_layer_sum(n_sp, nLayers, canopy_vol_frac, layer_id);

    Heightmidcrown_l[i] =
        Heightcrown_min_l[i] + (Height_max_l[i] - Heightcrown_min_l[i]) / 2.0;

    //determine the ratio between the mid height of the given species and the mid height of the layer.
    Heightmidcrown_r[i] = Heightmidcrown[i] / Heightmidcrown_l[i];

    // Calculate the sum of kL for all species in a layer
    kL_l[i] = k[i] * lai[i];
    kL_l[i] = f_get_layer_sum(n_sp, nLayers, kL_l, layer_id)[i];

    //Constant to partition light between species and to account for vertical canopy heterogeneity
    // (see Equations 2 and 3 of Forrester et al., 2014, Forest Ecosystems, 1:17)
    lambda_v[i] = 0.012306 +
        0.2366090 * k[i] * LAI[i] / kL_l[i] +
        0.029118 * Heightmidcrown_r[i] +
        0.608381 * k[i] * LAI[i] / kL_l[i] * Heightmidcrown_r[i];

    if (lai[i] == 0.0) {
      lambda_v[i] = 0.0;
    }
  }

// Normalize lambda_v to ensure the sum of all lambda_v values is 1
  lambdaV_l = f_get_layer_sum(n_sp, nLayers, lambda_v, layer_id);
  for (int i = 0; i < n_sp; i++) {
    if (lambdaV_l[i] != 0.0) {
      lambda_v[i] /= lambda_v[i] / lambdaV_l[i];
    }
  }

// Calculate the weighted kLS based on kL/sumkL
  for (int i = 0; i < n_sp; i++) {
    kLSweightedave[i] = k[i] * lai_sa_ratio[i] * k[i] * lai[i] / kL_l[i];
    kLSweightedave[i] =
        f_get_layer_sum(n_sp, nLayers, kLSweightedave, layer_id)[i];
    // the kLS should not be greater than 1 (based on the data used to fit the light model in Forrester et al. 2014)
    // This is because when there is a high k then LS is likely to be small
    if (kLSweightedave[i] > 1.0) {
      kLSweightedave[i] = 1.0;
    }

    // Constant to account for horizontal canopy heterogeneity such as gaps between trees and the change in zenith angle (and shading) with latitude and season (see Equations 2 and 5 of Forrester et al., 2014, Forest Ecosystems, 1:17)
    lambda_h[i] = 0.8285 +
        ((1.09498 - 0.781928 * kLSweightedave[i]) *
            pow(0.1, canopy_vol_frac[i])) -
        0.6714096 * pow(0.1, canopy_vol_frac[i]);

    if (solarAngle > 30.0) {
      lambda_h[i] += 0.00097 * pow(1.08259, solarAngle);
    }

    if (lai[i] == 0.0) {
      lambda_h[i] = 0.0;
    }
  }

// Calculate apar
  RADt = solar_rad * days_in_month;
  for (int i = 1; i <= nLayers; i++) {
    for (int j = 0; j < n_sp; j++) {
      if (layer_id[j] == i) {
        double apar_layer = RADt * (1.0 - exp(-kL_l[j]));
        apar[j] = apar_layer * lambda_h[j] * lambda_v[j];
        RADt -= apar_layer;
      }
    }
  }

// Calculate the proportion of above canopy apar absorbed by each species
  for (int i = 0; i < n_sp; i++) {
    fi[i] = apar[i] / (solar_rad * days_in_month);
  }

// Calculate the LAI above the given species for within canopy VPD calculations
  LAI_l = f_get_layer_sum(n_sp, nLayers, lambda_v, layer_id);

  // calculate the LAI of all layers above and part of the current layer if the species
  // is in the lower half of the layer then also take the proportion of the LAI above
  // the proportion is based on the Relative height of the mid crown
  for (int i = 0; i < n_sp; i++) {
    lai_above[i] = 0.0;

    for (int j = 0; j < n_sp; j++) {
      if (layer_id[j] < layer_id[i]) {
        lai_above[i] += lai[j];
      } else if (layer_id[j] == layer_id[i] &&
          Heightmidcrown_r[i] < 0.9999999999999) {
        lai_above[i] += lai[j] * (1.0 - Heightmidcrown_r[i]);
      }
    }
  }

  return [
    apar,
    lai_above,
    fi,
    lambda_v,
    lambda_h,
    canopy_vol_frac,
    layer_id,
    lai_sa_ratio
  ];
}

List<double> s_transpiration_3pgpjs(
  int n_sp,
  double solar_rad,
  double day_length,
  List<double> VPD_sp,
  List<double> BLcond,
  List<double> conduct_canopy,
  int days_in_month,
  List<double> Qa,
  List<double> Qb,
) {
  List<double> transp_veg = List<double>.filled(n_sp, 0.0); // Output

  List<double> netRad = List<double>.filled(n_sp, 0.0);
  List<double> defTerm = List<double>.filled(n_sp, 0.0);
  List<double> div = List<double>.filled(n_sp, 0.0);

  if (VPD_sp.reduce((sum, value) => sum + value) == 0.0) {
    transp_veg.fillRange(0, n_sp, 0.0);
  } else {
    for (int i = 0; i < n_sp; i++) {
      netRad[i] = Qa[i] + Qb[i] * (solar_rad * pow(10.0, 6.0) / day_length);
      defTerm[i] = rhoAir * lambda * (VPDconv * VPD_sp[i]) * BLcond[i];
      div[i] = conduct_canopy[i] * (1.0 + e20) + BLcond[i];
      transp_veg[i] = days_in_month *
          conduct_canopy[i] *
          (e20 * netRad[i] + defTerm[i]) /
          div[i] /
          lambda *
          day_length;
      transp_veg[i] = transp_veg[i] > 0.0 ? transp_veg[i] : 0.0;
    }
  }

  return transp_veg;
}

void s_transpiration_3pgmix(
  int n_sp,
  double solar_rad,
  double vpd_day,
  double day_length,
  int days_in_month,
  List<double> lai,
  List<double> fi,
  List<double> VPD_sp,
  List<double> aero_resist,
  List<double> conduct_canopy,
  double conduct_soil,
  List<double> Qa,
  List<double> Qb,
  List<double> transp_veg,
  double evapotra_soil,
) {
  List<double> netRad = List<double>.filled(n_sp, 0.0);
  List<double> defTerm = List<double>.filled(n_sp, 0.0);
  List<double> div = List<double>.filled(n_sp, 0.0);
  double lai_total;
  double netRad_so;
  double defTerm_so;
  double div_so;

  if (lai.reduce((sum, value) => sum + value) == 0.0) {
    transp_veg.fillRange(0, n_sp, 0.0);
  } else {
    for (int i = 0; i < n_sp; i++) {
      netRad[i] =
          (Qa[i] + Qb[i] * (solar_rad * pow(10.0, 6.0) / day_length)) * fi[i];
      defTerm[i] = rhoAir * lambda * (VPDconv * VPD_sp[i]) / aero_resist[i];
      div[i] = conduct_canopy[i] * (1.0 + e20) + 1.0 / aero_resist[i];
      transp_veg[i] = days_in_month *
          conduct_canopy[i] *
          (e20 * netRad[i] + defTerm[i]) /
          div[i] /
          lambda *
          day_length;
      if (lai[i] == 0.0) {
        transp_veg[i] = 0.0;
      }
    }
  }

  lai_total = lai.reduce((sum, value) => sum + value);

  if (lai_total > 0) {
    defTerm_so = rhoAir *
        lambda *
        (VPDconv * (vpd_day * exp(lai_total * (-ln2) / 5.0))) /
        (5.0 * lai_total);
    div_so = conduct_soil * (1.0 + e20) + 1.0 / (5.0 * lai_total);
  } else {
    defTerm_so =
        rhoAir * lambda * (VPDconv * (vpd_day * exp(lai_total * (-ln2) / 5.0)));
    div_so = conduct_soil * (1.0 + e20) + 1.0;
  }

  netRad_so = (Qa[0] + Qb[0] * (solar_rad * pow(10.0, 6.0) / day_length)) *
      (1.0 - fi.reduce((sum, value) => sum + value));
  evapotra_soil = days_in_month *
      conduct_soil *
      (e20 * netRad_so + defTerm_so) /
      div_so /
      lambda *
      day_length;
}

/// Returns dbh, basal_area, height, crown_length, crown_width, pFS, bias_scale in a List
List s_sizeDist_correct(
  int n_sp,
  List<double> age,
  List<double> stems_n,
  List<double> biom_tree,
  List<double> competition_total,
  List<double> lai,
  int correct_bias,
  int height_model,
  ParsS pars_s,
  ParsB pars_b,
  List<double> aWs,
  List<double> nWs,
  List<double> pfsPower,
  List<double> pfsConst,
  List<double> dbh, // Out also
  List<double> basal_area, // Out also
  List<double> height, // Out also
) {
  // Outputs
  List<double> crown_length = List<double>.filled(n_sp, 0.0);
  List<double> crown_width = List<double>.filled(n_sp, 0.0);
  List<double> pFS = List<double>.filled(n_sp, 0.0);
  List<List<double>> bias_scale =
      List.generate(15, (_) => List<double>.filled(n_sp, 0.0));

  double lai_total = 0;
  List<double> height_rel = List<double>.filled(n_sp, 0.0);

  // Additional variables for calculation distribution
  List<double> DWeibullScale = List<double>.filled(n_sp, 0.0);
  List<double> DWeibullShape = List<double>.filled(n_sp, 0.0);
  List<double> DWeibullLocation = List<double>.filled(n_sp, 0.0);
  List<double> wsWeibullScale = List<double>.filled(n_sp, 0.0);
  List<double> wsWeibullShape = List<double>.filled(n_sp, 0.0);
  List<double> wsWeibullLocation = List<double>.filled(n_sp, 0.0);
  List<double> Ex = List<double>.filled(n_sp, 0.0);
  List<double> Varx = List<double>.filled(n_sp, 0.0);
  List<double> CVdbhDistribution = List<double>.filled(n_sp, 0.0);
  List<double> CVwsDistribution = List<double>.filled(n_sp, 0.0);
  List<double> DrelBiaspFS = List<double>.filled(n_sp, 0.0);
  List<double> DrelBiasheight = List<double>.filled(n_sp, 0.0);
  List<double> DrelBiasBasArea = List<double>.filled(n_sp, 0.0);
  List<double> DrelBiasLCL = List<double>.filled(n_sp, 0.0);
  List<double> DrelBiasCrowndiameter = List<double>.filled(n_sp, 0.0);
  List<double> wsrelBias = List<double>.filled(n_sp, 0.0);

  List<double> dlocation = List<double>.filled(n_sp, 0.0);
  List<double> wslocation = List<double>.filled(n_sp, 0.0);
  List<double> DWeibullShape_gamma = List<double>.filled(n_sp, 0.0);
  List<double> wsWeibullShape_gamma = List<double>.filled(n_sp, 0.0);

  lai_total = lai.reduce((a, b) => a + b);

// Calculate the relative height
  double heightSum = 0.0;
  double stemsNSum = 0.0;

  for (int i = 0; i < n_sp; i++) {
    heightSum += height[i] * stems_n[i];
    stemsNSum += stems_n[i];
  }
  for (int i = 0; i < n_sp; i++) {
    height_rel[i] = height[i] / (heightSum / stemsNSum);
  }

// Check where all the locations are provided

  for (int i = 0; i < n_sp; i++) {
    if (pars_b.Dlocation0[i] == 0.0 &&
        pars_b.DlocationB[i] == 0.0 &&
        pars_b.Dlocationrh[i] == 0.0 &&
        pars_b.Dlocationt[i] == 0.0 &&
        pars_b.DlocationC[i] == 0.0) {
      dlocation[i] = 0.0;
    }
    if (pars_b.wslocation0[i] == 0.0 &&
        pars_b.wslocationB[i] == 0.0 &&
        pars_b.wslocationrh[i] == 0.0 &&
        pars_b.wslocationt[i] == 0.0 &&
        pars_b.wslocationC[i] == 0.0) {
      wslocation[i] = 0.0;
    }
  }

  if (correct_bias == 1) {
    // Calculate the DW scale
    for (int i = 0; i < n_sp; i++) {
      DWeibullScale[i] = exp(pars_b.Dscale0[i] +
          pars_b.DscaleB[i] * log(dbh[i]) +
          pars_b.Dscalerh[i] * log(height_rel[i]) +
          pars_b.Dscalet[i] * log(age[i]) +
          pars_b.DscaleC[i] * log(competition_total[i]));

      DWeibullShape[i] = exp(pars_b.Dshape0[i] +
          pars_b.DshapeB[i] * log(dbh[i]) +
          pars_b.Dshaperh[i] * log(height_rel[i]) +
          pars_b.Dshapet[i] * log(age[i]) +
          pars_b.DshapeC[i] * log(competition_total[i]));

      DWeibullShape_gamma[i] = f_gamma_dist(
          DWeibullShape.map((val) => 1.0 + 1.0 / val).toList(), n_sp)[i];

      DWeibullLocation[i] = exp(pars_b.Dlocation0[i] +
          pars_b.DlocationB[i] * log(dbh[i]) +
          pars_b.Dlocationrh[i] * log(height_rel[i]) +
          pars_b.Dlocationt[i] * log(age[i]) +
          pars_b.DlocationC[i] * log(competition_total[i]));

      if (dlocation[i] == 0.0) {
        DWeibullLocation[i] =
            (dbh[i] / 1.0 - 1.0) - DWeibullScale[i] * DWeibullShape_gamma[i];
      }

      if (DWeibullLocation[i] < 0.01) {
        DWeibullLocation[i] = 0.01;
      }

      Ex[i] = DWeibullLocation[i] + DWeibullScale[i] * DWeibullShape_gamma[i];

      // Now convert the Ex from Weibull scale to actual scale of diameter units in cm
      Varx[i] = pow(DWeibullScale[i], 2) *
          (f_gamma_dist(DWeibullShape.map((val) => 1.0 + 2.0 / val).toList(),
                  n_sp)[i] -
              pow(DWeibullShape_gamma[i], 2));

      CVdbhDistribution[i] = sqrt(Varx[i]) / Ex[i];

      for (int i = 0; i < n_sp; i++) {
        DrelBiaspFS[i] = 0.5 *
            pfsPower[i] *
            (pfsPower[i] - 1.0) *
            CVdbhDistribution[i] *
            CVdbhDistribution[i];
        DrelBiasheight[i] = 0.5 *
            pars_s.nHB[i] *
            (pars_s.nHB[i] - 1.0) *
            CVdbhDistribution[i] *
            CVdbhDistribution[i];
        DrelBiasBasArea[i] = 0.5 *
            2.0 *
            (2.0 - 1.0) *
            CVdbhDistribution[i] *
            CVdbhDistribution[i];
        DrelBiasLCL[i] = 0.5 *
            pars_s.nHLB[i] *
            (pars_s.nHLB[i] - 1.0) *
            CVdbhDistribution[i] *
            CVdbhDistribution[i];
        DrelBiasCrowndiameter[i] = 0.5 *
            pars_s.nKB[i] *
            (pars_s.nKB[i] - 1.0) *
            CVdbhDistribution[i] *
            CVdbhDistribution[i];

        DrelBiaspFS[i] = p_min_max(DrelBiaspFS, -0.5, 0.5)[i];
        DrelBiasheight[i] = p_min_max(DrelBiasheight, -0.5, 0.5)[i];
        DrelBiasBasArea[i] = p_min_max(DrelBiasBasArea, -0.5, 0.5)[i];
        DrelBiasLCL[i] = p_min_max(DrelBiasLCL, -0.5, 0.5)[i];
        DrelBiasCrowndiameter[i] =
            p_min_max(DrelBiasCrowndiameter, -0.5, 0.5)[i];
      }

      List<double> wsWeibullScale = List<double>.filled(n_sp, 0.0);
      List<double> wsWeibullShape = List<double>.filled(n_sp, 0.0);
      List<double> wsWeibullLocation = List<double>.filled(n_sp, 0.0);
      List<double> CVwsDistribution = List<double>.filled(n_sp, 0.0);

      wsWeibullScale[i] = exp(pars_b.wsscale0[i] +
          pars_b.wsscaleB[i] * log(dbh[i]) +
          pars_b.wsscalerh[i] * log(height_rel[i]) +
          pars_b.wsscalet[i] * log(age[i]) +
          pars_b.wsscaleC[i] * log(competition_total[i]));

      wsWeibullShape[i] = exp(pars_b.wsshape0[i] +
          pars_b.wsshapeB[i] * log(dbh[i]) +
          pars_b.wsshaperh[i] * log(height_rel[i]) +
          pars_b.wsshapet[i] * log(age[i]) +
          pars_b.wsshapeC[i] * log(competition_total[i]));

      double wsWeibullShape_gamma = f_gamma_dist(
          wsWeibullShape.map((val) => 1.0 + 1.0 / val).toList(), n_sp)[i];

      wsWeibullLocation[i] = exp(pars_b.wslocation0[i] +
          pars_b.wslocationB[i] * log(dbh[i]) +
          pars_b.wslocationrh[i] * log(height_rel[i]) +
          pars_b.wslocationt[i] * log(age[i]) +
          pars_b.wslocationC[i] * log(competition_total[i]));

      if (wslocation[i] == 0.0) {
        wsWeibullLocation[i] = (biom_tree[i].toInt() ~/ 10.0) -
            1.0 -
            wsWeibullScale[i] * wsWeibullShape_gamma;
      }

      if (wsWeibullLocation[i] < 0.01) {
        wsWeibullLocation[i] = 0.01;
      }

      Ex[i] = wsWeibullLocation[i] + wsWeibullScale[i] * wsWeibullShape_gamma;
      Varx[i] = wsWeibullScale[i] *
          wsWeibullScale[i] *
          (f_gamma_dist(wsWeibullShape.map((val) => 1.0 + 2.0 / val).toList(),
                  n_sp)[i] -
              wsWeibullShape_gamma * wsWeibullShape_gamma);
      CVwsDistribution[i] = sqrt(Varx[i]) / Ex[i];

      List<double> wsrelBias = List<double>.filled(n_sp, 0.0);

      wsrelBias[i] = 0.5 *
          (1.0 / nWs[i]) *
          (1.0 / nWs[i] - 1.0) *
          CVwsDistribution[i] *
          CVwsDistribution[i];
      wsrelBias[i] = p_min_max(wsrelBias, -0.5, 0.5)[i];
    }
  } else {
    DrelBiaspFS = List<double>.filled(n_sp, 0.0);
    DrelBiasBasArea = List<double>.filled(n_sp, 0.0);
    DrelBiasheight = List<double>.filled(n_sp, 0.0);
    DrelBiasLCL = List<double>.filled(n_sp, 0.0);
    DrelBiasCrowndiameter = List<double>.filled(n_sp, 0.0);
    wsrelBias = List<double>.filled(n_sp, 0.0);
  }

  // Correct for trees that have age 0 or are thinned (e.g. n_trees = 0)
  for (int i = 0; i < n_sp; i++) {
    if (age[i] == 0.0 || stems_n[i] == 0.0) {
      DrelBiaspFS[i] = 0.0;
      DrelBiasBasArea[i] = 0.0;
      DrelBiasheight[i] = 0.0;
      DrelBiasLCL[i] = 0.0;
      DrelBiasCrowndiameter[i] = 0.0;
      wsrelBias[i] = 0.0;
    }
  }

  // Correct for bias
  // Correct for bias
  for (int i = 0; i < n_sp; i++) {
    dbh[i] = pow(biom_tree[i] / aWs[i], 1.0 / nWs[i]) * (1.0 + wsrelBias[i]);
    basal_area[i] = pow(dbh[i], 2.0) /
        (4.0 * pi) *
        stems_n[i] /
        10000.0 *
        (1.0 + DrelBiasBasArea[i]);

    if (height_model == 1) {
      height[i] = pars_s.aH[i] *
          pow(dbh[i], pars_s.nHB[i]) *
          pow(competition_total[i], pars_s.nHC[i]) *
          (1.0 + DrelBiasheight[i]);
      crown_length[i] = pars_s.aHL[i] *
          pow(dbh[i], pars_s.nHLB[i]) *
          pow(lai_total, pars_s.nHLL[i]) *
          pow(competition_total[i], pars_s.nHLC[i]) *
          pow(height_rel[i], pars_s.nHLrh[i]) *
          (1.0 + DrelBiasLCL[i]);
    } else if (height_model == 2) {
      height[i] = 1.3 +
          pars_s.aH[i] * pow(exp(1.0), -pars_s.nHB[i] / dbh[i]) +
          pars_s.nHC[i] * competition_total[i] * dbh[i];
      crown_length[i] = 1.3 +
          pars_s.aHL[i] * pow(exp(1.0), -pars_s.nHLB[i] / dbh[i]) +
          pars_s.nHLC[i] * competition_total[i] * dbh[i];
    }

    crown_width[i] = pars_s.aK[i] *
        pow(dbh[i], pars_s.nKB[i]) *
        pow(height[i], pars_s.nKH[i]) *
        pow(competition_total[i], pars_s.nKC[i]) *
        pow(height_rel[i], pars_s.nKrh[i]) *
        (1.0 + DrelBiasCrowndiameter[i]);

    if (lai[i] == 0.0) {
      crown_width[i] = 0.0;
    }

    pFS[i] = pfsConst[i] * pow(dbh[i], pfsPower[i]) * (1.0 + DrelBiaspFS[i]);
  }

  // check that the height and LCL allometric equations have not predicted that height - LCL < 0
  // and if so reduce LCL so that height - LCL = 0 (assumes height allometry is more reliable than LCL
  for (int i = 0; i < n_sp; i++) {
    if (crown_length[i] > height[i]) {
      crown_length[i] = height[i];
    }
  }

  // Output the matrix of biases
  bias_scale[0] = DWeibullScale;
  bias_scale[1] = DWeibullShape;
  bias_scale[2] = DWeibullLocation;
  bias_scale[3] = wsWeibullScale;
  bias_scale[4] = wsWeibullShape;
  bias_scale[5] = wsWeibullLocation;
  bias_scale[6] = CVdbhDistribution;
  bias_scale[7] = CVwsDistribution;
  bias_scale[8] = wsrelBias;
  bias_scale[9] = DrelBiaspFS;
  bias_scale[10] = DrelBiasheight;
  bias_scale[11] = DrelBiasBasArea;
  bias_scale[12] = DrelBiasLCL;
  bias_scale[13] = DrelBiasCrowndiameter;
  bias_scale[14] = height_rel;

  return [dbh, basal_area, height, crown_length, crown_width, pFS, bias_scale];
}
