#pragma once
#include <cstdio>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "es/types.hpp"

namespace es {

/// Flat dotted-key configuration with physical unit suffixes.
///
/// Values may carry a unit, which is converted to SI on read:
///   length   nm um mm cm m          10um  -> 1e-5
///   voltage  V kV                   1.5kV -> 1500
///   current  pA nA uA mA A          80nA  -> 8e-8
///   pressure Pa kPa bar mbar        50Pa
///   flow     pL/s nL/s uL/s m3/s    2nL/s -> 2e-12
///   angle    deg rad                49.3deg
///   temp     K C                    25C   -> 298.15
/// A bare number is taken as already being in SI units.
///
/// Every key that is read is recorded, so unused_keys() reports typos -- a
/// silently ignored misspelled parameter is the most expensive kind of bug in a
/// parameter study.
class Config {
 public:
  void load(const std::string& path);
  void set(const std::string& key, const std::string& value);

  /// The non-"key=value" arguments (config file paths, --help, ...).  Call this
  /// first, load the files it names, and only THEN apply_cli -- command-line
  /// overrides must win over the file, not the other way round.
  static std::vector<std::string> positional_args(int argc, char** argv);
  /// Apply every "key=value" argument.
  void apply_cli(int argc, char** argv);

  bool has(const std::string& key) const;
  Real num(const std::string& key, Real fallback) const;
  int integer(const std::string& key, int fallback) const;
  bool flag(const std::string& key, bool fallback) const;
  std::string str(const std::string& key, const std::string& fallback) const;

  /// Keys that were never queried, excluding those under `ignore_prefixes`.
  /// Applications pass the sections they deliberately do not consume (es_field
  /// has no use for [beam], say), so that what remains really is a typo.
  std::vector<std::string> unused_keys(const std::vector<std::string>& ignore_prefixes = {}) const;
  void warn_about_unused(std::FILE* out,
                         const std::vector<std::string>& ignore_prefixes = {}) const;
  void dump(std::FILE* out) const;

  /// Parse one value with a unit suffix.  Exposed for testing.
  static Real parse_value(const std::string& text);

 private:
  std::map<std::string, std::string> kv_;
  mutable std::set<std::string> read_;
};

}  // namespace es
