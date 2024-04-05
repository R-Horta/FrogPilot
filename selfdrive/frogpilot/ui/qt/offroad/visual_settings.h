#pragma once

#include <set>

#include "selfdrive/ui/qt/offroad/settings.h"
#include "selfdrive/ui/ui.h"

class FrogPilotVisualsPanel : public FrogPilotListWidget {
  Q_OBJECT

public:
  explicit FrogPilotVisualsPanel(SettingsWindow *parent);

signals:
  void openParentToggle();

private:
  void hideToggles();
  void showEvent(QShowEvent *event);
  void updateCarToggles();
  void updateMetric();
  void updateState(const UIState &s);

  std::set<QString> alertVolumeControlKeys = {"DisengageVolume", "EngageVolume", "PromptDistractedVolume", "PromptVolume", "RefuseVolume", "WarningImmediateVolume", "WarningSoftVolume"};
  std::set<QString> customAlertsKeys = {};
  std::set<QString> customOnroadUIKeys = {"CustomPaths", "DeveloperUI", "LeadInfo"};
  std::set<QString> customThemeKeys = {"CustomColors", "CustomIcons", "CustomSignals", "CustomSounds"};
  std::set<QString> modelUIKeys = {"DynamicPathWidth", "HideLeadMarker", "LaneLinesWidth", "PathEdgeWidth", "PathWidth", "RoadEdgesWidth", "UnlimitedLength"};
  std::set<QString> qolKeys = {"CameraView", "DriverCamera"};
  std::set<QString> screenKeys = {};

  std::map<std::string, ParamControl*> toggles;

  Params params;

  bool hasAutoTune;
  bool hasBSM;
  bool hasOpenpilotLongitudinal;
  bool isMetric = params.getBool("IsMetric");
  bool isRelease;
  bool started;
};
