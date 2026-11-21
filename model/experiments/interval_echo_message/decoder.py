"""Terminal decision from one enrolled local outcome and a public lookup table."""


def decode(local_record, frozen_calibration):
    if type(local_record) is not int or not 0 <= local_record <= 12:
        raise ValueError("expected one local port or quiet")
    if (type(frozen_calibration) is not bytes or len(frozen_calibration) != 13
            or max(frozen_calibration) > 1):
        raise ValueError("expected thirteen binary choices")
    return frozen_calibration[local_record]
