class SizeDist {
  final String parameter;
  final Map<String, double> speciesValues;

  SizeDist({
    required this.parameter,
    required this.speciesValues,
  });
}

Map<String, SizeDist> prepareSizeDist(
    Map<String, dynamic>? sizeDist, List<String> spNames) {
  if (spNames.isEmpty) {
    throw Exception('spNames must be provided according to the species table.');
  }

  final Map<String, SizeDist> sizeDistOut = {};

  for (final spName in spNames) {
    sizeDistOut[spName] = SizeDist(
      parameter: '',
      speciesValues: {},
    );
  }

  final List<String> parameters =
      iSizeDist.map((item) => item['parameter'].toString()).toList();

  for (final spName in spNames) {
    sizeDistOut[spName] = SizeDist(
      parameter: '',
      speciesValues: {},
    );
    for (final parameter in parameters) {
      sizeDistOut[spName]!.speciesValues[parameter] = iSizeDist
          .firstWhere((item) => item['parameter'] == parameter)['default']
          .toDouble();
    }
  }

  if (sizeDist != null) {
    final List<String> sizeDistParameter = sizeDist.keys.toList();
    if (sizeDistParameter.isEmpty || sizeDistParameter.first != 'parameter') {
      throw Exception(
          'First column name of the size_dist table must correspond to: parameter');
    }

    final List<String> commonParameters =
        parameters.where((param) => sizeDistParameter.contains(param)).toList();

    if (commonParameters.length != sizeDistParameter.length - 1) {
      throw Exception(
          'size_dist input table must contain only parameters present in: ${parameters.join(',')}. Check `param_info` for more details.');
    }

    final List<String> spNamesReplace =
        spNames.where((spName) => sizeDist.containsKey(spName)).toList();

    for (final spName in spNamesReplace) {
      for (final parameter in commonParameters) {
        final value = sizeDist[parameter][spName];
        if (value != null) {
          sizeDistOut[spName]!.speciesValues[parameter] = value.toDouble();
        }
      }
    }
  }

  return sizeDistOut;
}
