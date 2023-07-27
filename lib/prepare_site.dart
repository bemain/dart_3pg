class Site {
  final double latitude;
  final double altitude;
  final int soilClass;
  final double asw_i;
  final double asw_min;
  final double asw_max;
  final String from;
  final String to;

  Site({
    required this.latitude,
    required this.altitude,
    required this.soilClass,
    required this.asw_i,
    required this.asw_min,
    required this.asw_max,
    required this.from,
    required this.to,
  });
}

Site prepareSite(Map<String, dynamic> site) {
  if (site.length != 8) {
    throw Exception('Site table shall contain exactly one row');
  }

  if (!site.containsKey('latitude') ||
      !site.containsKey('altitude') ||
      !site.containsKey('soil_class') ||
      !site.containsKey('asw_i') ||
      !site.containsKey('asw_min') ||
      !site.containsKey('asw_max') ||
      !site.containsKey('from') ||
      !site.containsKey('to')) {
    throw Exception(
        'Columns names of the site table must correspond to: latitude, altitude, soil_class, asw_i, asw_min, asw_max, from, to');
  }

  // Test for NA
  if (site.containsValue(null)) {
    throw Exception('Climate table should not contain null values');
  }

  // Prepare the time period
  final from = DateTime.parse('${site['from']}-01');
  final to = DateTime.parse('${site['to']}-01');

  if (from.isAfter(to)) {
    throw Exception('The start date is later than the end date');
  }

  final latitude = site['latitude'].toDouble();
  final altitude = site['altitude'].toDouble();
  final soilClass = site['soil_class'] as int;
  final asw_i = site['asw_i'].toDouble();
  final asw_min = site['asw_min'].toDouble();
  final asw_max = site['asw_max'].toDouble();
  final fromStr = site['from'].toString();
  final toStr = site['to'].toString();

  if (latitude > 90 || latitude < -90) {
    throw Exception('Latitude shall be within a range of [-90:90]');
  }

  if (altitude > 4000 || altitude < 0) {
    throw Exception('Altitude shall be within a range of [0:4000]');
  }

  if (soilClass < -1 || soilClass > 4) {
    throw Exception('Soil class shall be within a range of [-1:4]');
  }

  if (asw_i < 0) {
    throw Exception('ASW initial shall be greater than 0');
  }

  if (asw_min < 0) {
    throw Exception('ASW minimum shall be greater than 0');
  }

  if (asw_max < 0) {
    throw Exception('ASW maximum shall be greater than 0');
  }

  return Site(
    latitude: latitude,
    altitude: altitude,
    soilClass: soilClass,
    asw_i: asw_i,
    asw_min: asw_min,
    asw_max: asw_max,
    from: fromStr,
    to: toStr,
  );
}
