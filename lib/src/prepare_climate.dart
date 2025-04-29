import 'dart:math';

class Climate {
  final int year;
  final int month;
  final double tmpMin;
  final double tmpMax;
  final double tmpAve;
  final double prcp;
  final double srad;
  final double frostDays;
  final double vpdDay;
  final double co2;
  final double d13catm;

  Climate({
    required this.year,
    required this.month,
    required this.tmpMin,
    required this.tmpMax,
    required this.tmpAve,
    required this.prcp,
    required this.srad,
    required this.frostDays,
    required this.vpdDay,
    required this.co2,
    required this.d13catm,
  });
}

List<Climate> prepareClimate(
  List<Map<String, dynamic>> climate,
  String from,
  String to,
) {
  if (!climate.every((item) =>
      item.containsKey('tmp_min') &&
      item.containsKey('tmp_max') &&
      item.containsKey('prcp') &&
      item.containsKey('srad') &&
      item.containsKey('frost_days'))) {
    throw Exception(
        'Climate table must include the following columns: tmp_min, tmp_max, prcp, srad, frost_days');
  }

  if (climate.any((item) => item.containsValue(null))) {
    throw Exception('Climate table should not contain null values');
  }

  final fromDate = DateTime.parse(from + '-01');
  final toDate = DateTime.parse(to + '-01');

  if (fromDate.isAfter(toDate)) {
    throw Exception('The start date is later than the end date');
  }

  List<Climate> climateData = [];

  if (climate.length == 12) {
    final int nYears = toDate.year - fromDate.year + 1;
    final int monthI = fromDate.month;
    final int monthE = toDate.month;

    for (int i = 0; i < nYears; i++) {
      for (int j = 0; j < 12; j++) {
        climateData.add(Climate(
          year: fromDate.year + i,
          month: j + 1,
          tmpMin: climate[j]['tmp_min'],
          tmpMax: climate[j]['tmp_max'],
          tmpAve: climate[j]['tmp_ave'] ??
              (climate[j]['tmp_min'] + climate[j]['tmp_max']) / 2,
          prcp: climate[j]['prcp'],
          srad: climate[j]['srad'],
          frostDays: climate[j]['frost_days'],
          vpdDay: getVpd(climate[j]['tmp_min'], climate[j]['tmp_max']),
          co2: climate[j]['co2'] ?? 350,
          d13catm: climate[j]['d13catm'] ?? -7.1,
        ));
      }
    }

    if (monthI > 1) {
      climateData = climateData.skip(monthI - 1).toList();
    }

    if (monthE < 12) {
      climateData =
          climateData.take(climateData.length - (12 - monthE)).toList();
    }
  } else {
    if (!climate.every(
        (item) => item.containsKey('year') && item.containsKey('month'))) {
      throw Exception(
          'Climate table must include year and month for subsetting.');
    }

    climateData = climate
        .map((item) => Climate(
              year: item['year'],
              month: item['month'],
              tmpMin: item['tmp_min'],
              tmpMax: item['tmp_max'],
              tmpAve:
                  item['tmp_ave'] ?? (item['tmp_min'] + item['tmp_max']) / 2,
              prcp: item['prcp'],
              srad: item['srad'],
              frostDays: item['frost_days'],
              vpdDay: getVpd(item['tmp_min'], item['tmp_max']),
              co2: item['co2'] ?? 350,
              d13catm: item['d13catm'] ?? -7.1,
            ))
        .toList();

    final List<DateTime> dates =
        climateData.map((item) => DateTime(item.year, item.month)).toList();

    if (fromDate.isBefore(dates.first) || toDate.isAfter(dates.last)) {
      throw Exception(
          'Requested time period is outside of provided dates in climate table.');
    }

    climateData = climateData
        .where((item) =>
            item.year >= fromDate.year &&
            item.year <= toDate.year &&
            item.month >= fromDate.month &&
            item.month <= toDate.month)
        .toList();
  }

  return climateData;
}

double getVpd(double tMin, double tMax) {
  final vpdMin = 6.10780 * exp(17.2690 * tMin / (237.30 + tMin));
  final vpdMax = 6.10780 * exp(17.2690 * tMax / (237.30 + tMax));

  return (vpdMax - vpdMin) / 2;
}
