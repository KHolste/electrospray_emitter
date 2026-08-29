#include "es/feed.hpp"

#include <cmath>
#include <limits>

#include "es/constants.hpp"

namespace es {

using constants::pi;

const char* to_string(PressureMode m) {
  switch (m) {
    case PressureMode::Direct: return "Direct";
    case PressureMode::Budget: return "Budget";
  }
  return "?";
}

const char* to_string(FeedStatus s) {
  switch (s) {
    case FeedStatus::NotAttempted: return "NotAttempted";
    case FeedStatus::Ok: return "Ok";
    case FeedStatus::MissingFeedInput: return "MissingFeedInput";
    case FeedStatus::ChannelGeometryInvalid: return "ChannelGeometryInvalid";
    case FeedStatus::MissingLiquidProperty: return "MissingLiquidProperty";
    case FeedStatus::NotLaminar: return "NotLaminar";
    case FeedStatus::EntranceLengthNotShort: return "EntranceLengthNotShort";
  }
  return "?";
}

const char* explain(FeedStatus s) {
  switch (s) {
    case FeedStatus::NotAttempted:
      return "Der Druckhaushalt wurde nicht ausgewertet.";
    case FeedStatus::Ok:
      return "Der Druckhaushalt ist ausgewertet und die Gueltigkeitsgrenzen der "
             "Poiseuille-Beziehung sind eingehalten.";
    case FeedStatus::MissingFeedInput:
      return "Eine Pflichtangabe des Druckhaushalts fehlt oder ist nicht physikalisch.  Es "
             "wird KEIN Ersatzwert gesetzt: der Reservoirdruck und der Volumenstrom sind "
             "Eingaben und werden hier nicht berechnet.";
    case FeedStatus::ChannelGeometryInvalid:
      return "Radius oder Laenge des Zulaufkanals sind nicht positiv.";
    case FeedStatus::MissingLiquidProperty:
      return "Der Stoffdatensatz traegt keine Viskositaet oder keine Dichte.  Ohne sie gibt "
             "es weder einen viskosen noch einen hydrostatischen Term, und ein Nullwert "
             "waere eine Behauptung.";
    case FeedStatus::NotLaminar:
      return "Die Reynoldszahl uebersteigt die laminare Schranke.  Die Hagen-Poiseuille-"
             "Beziehung gilt dort nicht; das Ergebnis wird deshalb nicht als gueltig "
             "ausgewiesen.";
    case FeedStatus::EntranceLengthNotShort:
      return "Die hydrodynamische Einlauflaenge ist nicht klein gegen die Kanallaenge.  Die "
             "Stroemung ist dann nicht ueber die ganze Laenge ausgebildet und der "
             "geschlossene Druckverlust ist nicht der ganze Verlust.";
  }
  return "?";
}

// ---------------------------------------------------------------------------

Real FeedChannel::hydraulic_resistance(Real mu) const {
  if (!(radius > 0.0) || !(length > 0.0)) return std::numeric_limits<Real>::quiet_NaN();
  const Real r2 = radius * radius;
  return 8.0 * mu * length / (pi * r2 * r2);
}

Real FeedChannel::mean_velocity(Real Q) const {
  if (!(radius > 0.0)) return std::numeric_limits<Real>::quiet_NaN();
  return Q / (pi * radius * radius);
}

Real FeedChannel::centreline_velocity(Real Q) const { return 2.0 * mean_velocity(Q); }

Real FeedChannel::velocity_at(Real r, Real Q) const {
  if (!(radius > 0.0)) return std::numeric_limits<Real>::quiet_NaN();
  const Real x = r / radius;
  return centreline_velocity(Q) * (1.0 - x * x);
}

Real FeedChannel::wall_shear_stress(Real mu, Real Q) const {
  if (!(radius > 0.0)) return std::numeric_limits<Real>::quiet_NaN();
  return 4.0 * mu * mean_velocity(Q) / radius;
}

Real FeedChannel::reynolds(Real rho, Real mu, Real Q) const {
  if (!(radius > 0.0) || !(mu > 0.0)) return std::numeric_limits<Real>::quiet_NaN();
  return rho * std::abs(mean_velocity(Q)) * 2.0 * radius / mu;
}

Real FeedChannel::entrance_length(Real rho, Real mu, Real Q) const {
  return 0.06 * reynolds(rho, mu, Q) * 2.0 * radius;
}

// ---------------------------------------------------------------------------

PressureBudget solve_pressure_budget(const FeedRequest& q, const LiquidProperties& liquid,
                                     Real contact_radius) {
  PressureBudget b;
  b.mode = q.mode;
  b.gamma_over_a = (contact_radius > 0.0) ? liquid.gamma / contact_radius : 0.0;

  auto finish = [&]() {
    b.Pi = (b.gamma_over_a > 0.0) ? b.delta_p_exit / b.gamma_over_a : 0.0;
    b.within_capillary_range = (b.gamma_over_a > 0.0) &&
                               (std::abs(b.delta_p_exit) <= 2.0 * b.gamma_over_a);
    return b;
  };

  if (q.contact_angle_requested) {
    b.status = FeedStatus::MissingFeedInput;
    b.message =
        "Ein Kontaktwinkel wurde gesetzt.  Der Zulaufkanal ist in diesem Modell VOLL "
        "gefuellt; dort existiert keine zweite freie Oberflaeche, deren Young-Winkel "
        "eingesetzt werden koennte.  Die einzige freie Oberflaeche ist der an der "
        "Austrittskante gepinnte Meniskus, und dessen Kapillardruck rechnet P3a/P3b -- "
        "er darf hier nicht ein zweites Mal addiert werden.";
    b.delta_p_exit = std::numeric_limits<Real>::quiet_NaN();
    return finish();
  }

  if (q.mode == PressureMode::Direct) {
    b.delta_p_exit = q.delta_p_exit_direct;
    b.driving = std::numeric_limits<Real>::quiet_NaN();
    b.hydrostatic = std::numeric_limits<Real>::quiet_NaN();
    b.viscous = std::numeric_limits<Real>::quiet_NaN();
    b.status = FeedStatus::Ok;
    b.message = "Direkte Vorgabe von delta_p_exit; kein Haushalt gerechnet.  Die einzelnen "
                "Terme sind deshalb nicht bekannt und stehen als nan.";
    return finish();
  }

  // --- the budget -----------------------------------------------------------
  const Real mu = liquid.documented_only.mu;
  const Real rho = liquid.rho;
  if (!(q.channel.radius > 0.0) || !(q.channel.length > 0.0)) {
    b.status = FeedStatus::ChannelGeometryInvalid;
    b.message = explain(b.status);
    b.delta_p_exit = std::numeric_limits<Real>::quiet_NaN();
    return finish();
  }
  const bool needs_mu = (q.Q != 0.0);
  const bool needs_rho = (q.gravity_axial != 0.0) && (q.z_exit != q.z_reservoir);
  if ((needs_mu && !(mu > 0.0)) || (needs_rho && !(rho > 0.0))) {
    b.status = FeedStatus::MissingLiquidProperty;
    b.message = explain(b.status);
    b.delta_p_exit = std::numeric_limits<Real>::quiet_NaN();
    return finish();
  }

  b.driving = q.p_reservoir - q.p_vacuum;
  // p(z) = p_res + rho g_z (z - z_res), so the loss on the way up is
  // p_res - p(z_exit) = -rho g_z (z_exit - z_res).
  b.hydrostatic = -rho * q.gravity_axial * (q.z_exit - q.z_reservoir);
  b.hydraulic_resistance = q.channel.hydraulic_resistance(mu);
  b.viscous = b.hydraulic_resistance * q.Q;
  b.delta_p_exit = b.driving - b.hydrostatic - b.viscous;

  b.mean_velocity = q.channel.mean_velocity(q.Q);
  b.wall_shear_stress = q.channel.wall_shear_stress(mu, q.Q);
  b.reynolds = (mu > 0.0) ? q.channel.reynolds(rho, mu, q.Q) : 0.0;
  b.entrance_length = (mu > 0.0) ? q.channel.entrance_length(rho, mu, q.Q) : 0.0;
  b.entrance_fraction = b.entrance_length / q.channel.length;

  if (!(b.reynolds < feed::kReynoldsLaminar)) {
    b.status = FeedStatus::NotLaminar;
    b.message = explain(b.status);
    return finish();
  }
  if (!(b.entrance_fraction < feed::kEntranceFraction)) {
    b.status = FeedStatus::EntranceLengthNotShort;
    b.message = explain(b.status);
    return finish();
  }
  b.status = FeedStatus::Ok;
  b.message = explain(b.status);
  return finish();
}

void PressureBudget::print(std::FILE* out) const {
  std::fprintf(out, "  Druckhaushalt am Austritt (%s): %s\n", to_string(mode), to_string(status));
  std::fprintf(out, "    delta_p_exit          = %.6g Pa   (Pi = %.4f)\n", delta_p_exit, Pi);
  if (mode == PressureMode::Budget) {
    std::fprintf(out, "      p_res - p_vak       = %+.6g Pa\n", driving);
    std::fprintf(out, "      - hydrostatisch     = %+.6g Pa\n", -hydrostatic);
    std::fprintf(out, "      - viskos            = %+.6g Pa\n", -viscous);
    std::fprintf(out, "    Widerstand R_h        = %.6g Pa s/m^3\n", hydraulic_resistance);
    std::fprintf(out, "    mittlere Geschwindigkeit %.6g m/s, Wandschubspannung %.6g Pa\n",
                 mean_velocity, wall_shear_stress);
    std::fprintf(out, "    Re = %.4g (Grenze %.0f), Einlauflaenge %.4g m = %.4g L\n", reynolds,
                 feed::kReynoldsLaminar, entrance_length, entrance_fraction);
  }
  std::fprintf(out, "    gamma/a = %.6g Pa, im Kapillarbereich |dp| <= 2 gamma/a: %s\n",
               gamma_over_a, within_capillary_range ? "ja" : "NEIN");
  if (!message.empty()) std::fprintf(out, "    %s\n", message.c_str());
}

}  // namespace es
