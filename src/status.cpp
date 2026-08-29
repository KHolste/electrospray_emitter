#include "es/status.hpp"

namespace es {

const char* to_string(SolveStatus s) {
  switch (s) {
    case SolveStatus::Converged: return "converged";
    case SolveStatus::NotConverged: return "not_converged";
    case SolveStatus::VoltageNotBracketed: return "voltage_not_bracketed";
    case SolveStatus::VoltageMismatch: return "voltage_mismatch";
    case SolveStatus::ShapeIntegrationFailed: return "shape_integration_failed";
    case SolveStatus::NoStaticFoldFound: return "no_static_fold_found";
    default: return "not_attempted";
  }
}

const char* explain(SolveStatus s) {
  switch (s) {
    case SolveStatus::Converged:
      return "Loesung konvergiert und trifft die angeforderte Vorgabe.";
    case SolveStatus::NotConverged:
      return "Die Form-Feld-Iteration hat ihre Toleranz nicht erreicht.";
    case SolveStatus::VoltageNotBracketed:
      return "Die angeforderte Spannung liegt ausserhalb des erreichbaren Astes. "
             "Unterhalb der kleinsten Astspannung existiert bei diesem Speisedruck "
             "kein Gleichgewicht, oberhalb der Faltenspannung keine statische Loesung.";
    case SolveStatus::VoltageMismatch:
      return "Die Suche endete bei einer anderen Spannung als angefordert. "
             "Die gelieferte Form gehoert NICHT zur angeforderten Spannung.";
    case SolveStatus::ShapeIntegrationFailed:
      return "Die Young-Laplace-Integration hat die Kontaktlinie nicht erreicht "
             "(Oberflaeche wird senkrecht, r nicht mehr monoton).";
    case SolveStatus::NoStaticFoldFound:
      return "Der Ast besitzt keinen inneren Umkehrpunkt, unter dem gesucht werden koennte.";
    default:
      return "Es wurde nichts geloest.";
  }
}

const char* to_string(FoldStatus s) {
  switch (s) {
    case FoldStatus::Found: return "found";
    case FoldStatus::TooFewPoints: return "too_few_points";
    case FoldStatus::Monotone: return "monotone";
    case FoldStatus::MaximumAtBoundary: return "maximum_at_boundary";
    default: return "not_attempted";
  }
}

const char* explain(FoldStatus s) {
  switch (s) {
    case FoldStatus::Found:
      return "Innerer Umkehrpunkt des statischen Astes nachgewiesen.";
    case FoldStatus::TooFewPoints:
      return "Weniger als drei konvergierte Astpunkte: ein Umkehrpunkt ist damit "
             "nicht feststellbar.";
    case FoldStatus::Monotone:
      return "Die Spannung ist entlang des gesamten Astes monoton. Es gibt keinen "
             "Umkehrpunkt; der Ast wurde vermutlich nicht weit genug verfolgt.";
    case FoldStatus::MaximumAtBoundary:
      return "Das Maximum liegt am Rand des verfolgten Bereichs und ist damit kein "
             "nachgewiesener innerer Umkehrpunkt. Bereich erweitern.";
    default:
      return "Es wurde kein Ast ausgewertet.";
  }
}

}  // namespace es
