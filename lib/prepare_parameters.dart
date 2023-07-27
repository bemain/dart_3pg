class Parameters {
  final String parameter;
  final Map<String, dynamic> speciesValues;

  Parameters({
    required this.parameter,
    required this.speciesValues,
  });
}

Map<String, Parameters> prepareParameters(
    Map<String, dynamic>? parameters, List<String> spNames) {
  if (spNames.isEmpty) {
    throw Exception('spNames must be provided according to the species table.');
  }

  final Map<String, Parameters> parametersOut = {};

  for (final spName in spNames) {
    parametersOut[spName] = Parameters(
      parameter: '',
      speciesValues: {},
    );
  }

  final List<String> parameterList =
      iParameters.map((item) => item['parameter'].toString()).toList();

  for (final spName in spNames) {
    parametersOut[spName] = Parameters(
      parameter: '',
      speciesValues: {},
    );
    for (final parameter in parameterList) {
      parametersOut[spName]!.speciesValues[parameter] = iParameters
          .firstWhere((item) => item['parameter'] == parameter)['default'];
    }
  }

  if (parameters != null) {
    final List<String> parameterColumns = parameters.keys.toList();
    if (parameterColumns.isEmpty || parameterColumns.first != 'parameter') {
      throw Exception(
          'First column name of the parameters table must correspond to: parameter');
    }

    final List<String> commonParameters = parameterList
        .where((param) => parameterColumns.contains(param))
        .toList();

    if (commonParameters.length != parameterColumns.length - 1) {
      throw Exception(
          'Parameter input table must contain only parameters present in: ${parameterList.join(',')}. Check `param_info` for more details.');
    }

    final List<String> spNamesReplace =
        spNames.where((spName) => parameters.containsKey(spName)).toList();

    for (final spName in spNamesReplace) {
      for (final parameter in commonParameters) {
        final value = parameters[parameter][spName];
        if (value != null) {
          parametersOut[spName]!.speciesValues[parameter] = value;
        }
      }
    }
  }

  return parametersOut;
}
