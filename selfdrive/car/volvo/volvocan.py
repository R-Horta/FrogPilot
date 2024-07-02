from openpilot.selfdrive.car.volvo.values import Buttons
from openpilot.common.conversions import Conversions as CV

def create_button_msg(packer, resume=False, cancel=False, bus=0):
  # TODO: validate
  msg = {
    "ACCOnOffBtn": cancel,
    "ACCOnOffBtnInv": not cancel,
    "ACCResumeBtn": resume,
    "ACCResumeBtnInv": not resume,
  }
  return packer.make_can_msg("CCButtons", bus, msg)


def create_acc_state_msg(packer):
  msg = {
    "ACC_Check": 1,
  }
  return packer.make_can_msg("FSM3", 0, msg)


def create_lkas_state_msg(packer, steering_angle: float, stock_values: dict):
  # Manipulate data from servo to FSM
  # Set LKATorque and LKAActive to zero otherwise LKA will be disabled. (Check dbc)
  msg = {
    "LKATorque": 0,
    "SteeringAngleServo": steering_angle,
    "byte0": stock_values["byte0"],
    "byte4": stock_values["byte4"],
    "byte7": stock_values["byte7"],
    "LKAActive": int(stock_values["LKAActive"]) & 0xF5,
    "SteeringWheelRateOfChange": stock_values["SteeringWheelRateOfChange"],
  }
  return packer.make_can_msg("PSCM1", 2, msg)


def calculate_lka_checksum(dat: bytes) -> int:
  # Input: dat byte array, and fingerprint
  # Steering direction = 0 -> 3
  # TrqLim = 0 -> 255
  # Steering angle request = -360 -> 360

  # Extract LKAAngleRequest, LKADirection and Unknown
  steer_angle_request = ((dat[3] & 0x3F) << 8) + dat[4]
  steering_direction_request = dat[5] & 0x03
  trqlim = dat[2]

  # Sum of all bytes, carry ignored.
  s = (trqlim + steering_direction_request + steer_angle_request + (steer_angle_request >> 8)) & 0xFF
  # Checksum is inverted sum of all bytes
  return s ^ 0xFF


def create_lka_msg(packer, apply_steer: float, steer_direction: int):
  values = {
    "LKAAngleReq": apply_steer,
    "LKASteerDirection": steer_direction,
    "TrqLim": 0,

    # car specific parameters
    "SET_X_22": 0x25, # Test these values: 0x24, 0x22
    "SET_X_02": 0,    # Test 0x00, 0x02
    "SET_X_10": 0x10, # Test 0x10, 0x1c, 0x18, 0x00
    "SET_X_A4": 0xa7, # Test 0xa4, 0xa6, 0xa5, 0xe5, 0xe7
  }

  # calculate checksum
  dat = packer.make_can_msg("FSM2", 0, values)[2]
  values["Checksum"] = calculate_lka_checksum(dat)

  return packer.make_can_msg("FSM2", 0, values)

def create_button_cmd(packer, CP, counter, button):

  can = int(button == Buttons.CANCEL)
  res = int(button == Buttons.RESUME)
  inc = int(button == Buttons.SET_PLUS)
  dec = int(button == Buttons.SET_MINUS)

  values = {
    "ACCOnOffBtn": can,
    "ACCOnOffBtnInv": (can + 1) % 2,

    "ACCSetBtn": inc,
    "ACCSetBtnInv": (inc + 1) % 2,

    "ACCResumeBtn": res,
    "ACCResumeBtnInv": (res + 1) % 2,

    "ACCMinusBtn": dec,
    "ACCMinusBtnInv": (dec + 1) % 2,

    "TimeGapDecreaseBtn": 0,
    "TimeGapDecreaseBtnInv": 1,

    "TimeGapIncreaseBtn": 0,
    "TimeGapIncreaseBtnInv": 1,
    }

  return packer.make_can_msg("CCButtons", 0, values)

def create_volvo_acc_spam_command(packer, controller, CS, slcSet, Vego, frogpilot_toggles, accel):
  cruiseBtn = Buttons.NONE

  MS_CONVERT = CV.MS_TO_KPH if frogpilot_toggles.is_metric else CV.MS_TO_MPH

  speedSetPoint = int(round(CS.out.cruiseState.speed * MS_CONVERT))
  slcSet = int(round(slcSet * MS_CONVERT))

  #if not frogpilot_toggles.experimentalMode:
  if slcSet + 5 < Vego * MS_CONVERT:
    slcSet = slcSet - 10 # 10 lower to increase deceleration until with 5
  #else:
  #  slcSet = int(round((Vego + 5 * accel) * MS_CONVERT))
  
  if frogpilot_toggles.is_metric: # Default is by 5 kph
    slcSet = int(round(slcSet/5.0)*5.0)
    speedSetPoint = int(round(speedSetPoint/5.0)*5.0)

  if slcSet < speedSetPoint and speedSetPoint > (32 if frogpilot_toggles.is_metric else 20):
    cruiseBtn = Buttons.SET_MINUS
  elif slcSet > speedSetPoint:
    cruiseBtn = Buttons.SET_PLUS
  else:
    cruiseBtn = Buttons.NONE

  if (cruiseBtn != Buttons.NONE):
    return [create_button_cmd(packer, controller.CP, controller.frame // 10, cruiseBtn)]
  else:
    return []
