import 'package:collection/collection.dart';
import 'prepare_climate.dart';
import 'prepare_parameters.dart';
import 'prepare_site.dart';
import 'prepare_sizeDist.dart';
import 'prepare_thinning.dart';

Map<String, dynamic> prepareInput({
  required Map<String, dynamic> site,
  required List<Map<String, dynamic>> species,
  List<Map<String, dynamic>>? climate,
  List<Map<String, dynamic>>? thinning,
  List<Map<String, dynamic>>? parameters,
  List<Map<String, dynamic>>? sizeDist,
  Map<String, dynamic>? settings,
}) {
  // Site
  Site preparedSite = prepareSite(site: site);

  // Species
  List<Map<String, dynamic>> preparedSpecies = prepareSpecies(species: species);

  // Settings
  Map<String, dynamic> defaultSettings = {
    'light_model': 1,
    'transp_model': 1,
    'phys_model': 1,
    'height_model': 1,
    'correct_bias': 0,
    'calculate_d13c': 0
  };
  if (settings != null) {
    defaultSettings.addAll(settings);
  }

  // Climate
  if (defaultSettings['calculate_d13c'] == 1) {
    if (!(climate![0].containsKey('co2') &&
        climate[0].containsKey('d13catm'))) {
      throw Exception(
          'Please provide forcing data for co2 and d13catm in climate, if calculate_d13c = 1');
    }
  }
  List<Map<String, dynamic>> preparedClimate = prepareClimate(
      climate: climate!, from: preparedSite['from'], to: preparedSite['to']);

  // Thinning
  List<List<List<double>>> preparedThinning = prepareThinning(
    thinning: thinning,
    spNames: preparedSpecies.map((row) => row['species'] as String).toList(),
  );

  // Parameters
  List<Map<String, dynamic>> preparedParameters = prepareParameters(
    parameters: parameters,
    spNames: preparedSpecies.map((row) => row['species'] as String).toList(),
  );

  // Size distribution
  if (defaultSettings['correct_bias'] == 1 && sizeDist == null) {
    throw Exception(
        'Please provide size_dist table or change the setting to size_dist = 0');
  }
  List<Map<String, dynamic>> preparedSizeDist = prepareSizeDist(
    sizeDist: sizeDist,
    spNames: preparedSpecies.map((row) => row['species'] as String).toList(),
  );

  // return the checked output
  Map<String, dynamic> output = {
    'site': preparedSite,
    'species': preparedSpecies,
    'climate': preparedClimate,
    'thinning': preparedThinning,
    'parameters': preparedParameters,
    'size_dist': preparedSizeDist,
    'settings': defaultSettings,
  };

  return output;
}
