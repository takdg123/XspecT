from astropy.time.formats import TimeFromEpoch
from astropy.time import Time
import datetime


class TimeFermiSec(TimeFromEpoch):
    """Represents the number of seconds elapsed since Jan 1, 2001 00:00:00 UTC including leaps seconds"""

    @property
    def value(self):
        return super(TimeFermiSec, self).value()

    name = 'fermi'
    unit = 1.0 / 86400  # in days (1 day == 86400 seconds)
    epoch_val = '2001-01-01 00:01:04.184'
    epoch_val2 = None
    epoch_scale = 'tt'  # Scale for epoch_val class attribute
    epoch_format = 'iso'  # Format for epoch_val class attribute


def round_half_to_nearest_even(num):
    """Round the given number to the nearest even integer value.

    Parameters
    ==========
    :param num: floating point number
    """
    n = int(num)
    v = abs(num - n)
    if v > 0.5 or (v == 0.5 and n % 2):
        return n + 1 - 2 * int(n < 0)
    else:
        return n


def hms_to_fraction_of_day(value):
    """The fraction of day as computed by the original pipeline code.

    Parameters
    ==========
    :param value: datetime object"""
    result = round_half_to_nearest_even(((value.hour * 3600 + value.minute * 60 + value.second) / 86400) * 1000)
    return min(result, 999)


def fraction_of_day_to_hms(value):
    """Returns the hour, minute, second for a given fraction of day"""
    s = int((value / 1000) * 86400)
    h = s // 3600
    s -= h * 3600
    m = s // 60
    s -= m * 60
    return h, m, s


class Met(object):
    """Represents the Fermi timescale"""

    # Mission Elapsed Time (Number of seconds since 2001-01-01 00:00:00 UTC)

    def __init__(self, secs):
        """ Creates a Met object with the time set to the number of seconds since Jan 1, 2001 00:00:00 UTC including the
         leap seconds"""
        if secs < 0:
            raise Exception("Time before GBM mission epoch")
        self.__time = Time(secs, format='fermi')

    @property
    def met(self):
        return self.__time.fermi

    # Astropy Time

    @property
    def time(self):
        """Returns the value as a astropy.Time object"""
        return self.__time

    @classmethod
    def from_time(cls, atime):
        """Creates a new Met object set to the given astropy.Time object"""
        obj = cls(0)
        obj.__time = atime
        if obj.met < 0:
            raise Exception("Time before GBM mission epoch")
        return obj

    # Python's datetime

    @property
    def datetime(self):
        """Returns a datetime object set to the correct UTC value."""
        try:
            return self.__time.utc.to_datetime(datetime.timezone.utc)
        except ValueError:
            # Repeat last met for a leap second
            return Met(self.met-1).datetime

    @classmethod
    def from_datetime(cls, dt):
        """Create a Met object with its value set to the given datetime.datetime value"""
        return cls.from_time(Time(dt, format='datetime'))

    # Unix timestamp (Number of seconds since 1970-01-01 00:00:00 UTC (leap seconds are ignored))

    @property
    def unix(self):
        """Returns the value as the number of seconds since Jan 1, 1970 00:00:00 with the leap seconds removed."""
        return self.datetime.timestamp()

    @classmethod
    def from_unix(cls, unix):
        """Creates a Met object with its value set to the given Unix time value"""
        return cls.from_datetime(datetime.datetime.utcfromtimestamp(unix))

    # GPS timestamp (Number of seconds since Jan 6, 1980 00:00:00 UTC (leap seconds are ignored))

    @property
    def gps(self):
        """Returns the value as the number of seconds since Jan 6, 1980 00:00:00 (leap seconds are ignored)"""
        return self.__time.gps

    @classmethod
    def from_gps(cls, gps):
        """Creates a Met object with its value set to the given GPS time value"""
        return cls.from_time(Time(gps, format='gps'))

    # Julian date

    @property
    def jd(self):
        """Returns the value as the Julian Date"""
        return self.__time.jd

    @classmethod
    def from_jd(cls, jd):
        """Creates new Met object with the value of the given Julian Date"""
        return cls.from_time(Time(jd, format='jd'))

    # Modified Julian Date

    @property
    def mjd(self):
        """Returns the value as a Modified Julian Date"""
        return self.__time.mjd

    @classmethod
    def from_mjd(cls, mjd):
        return cls.from_time(Time(mjd, format='mjd'))

    # GBM Burst Number (YYMMDDFFF)

    @property
    def bn(self):
        """Returns the MET value as a burst number string in the form of YYMMDDFFF."""

        # Adjust to match a known bug in the old pipeline software
        adj_met = self.met
        if 157766399.0 < adj_met < 252460801.0:
            adj_met += 1
        elif 252460800.0 < adj_met <= 253497600.0:
            adj_met += 2

        # To ensure compatibility with the number produced by the pipeline, we are doing it the inefficient way
        m = Met(adj_met)
        utc_val = m.datetime
        fraction = hms_to_fraction_of_day(utc_val)

        return "{}{:03d}".format(utc_val.strftime("%y%m%d"), fraction)

    @classmethod
    def from_bn(cls, bn):
        """Create a Met object with its value set to the equivalent of YYMMDDFFF string"""
        dt = datetime.datetime.strptime(bn[:6], '%y%m%d')
        hms = fraction_of_day_to_hms(int(bn[6:]))
        dt = datetime.datetime(dt.year, dt.month, dt.day, hms[0], hms[1], hms[2], tzinfo=datetime.timezone.utc)
        obj = cls.from_datetime(dt)

        # Adjust to match a known bug in the old pipeline software
        adj_met = obj.met
        if 157766400.0 < adj_met < 252460802.0:
            adj_met -= 1
        elif 252460802.0 < adj_met <= 253497602.0:
            adj_met -= 2

        return Met(adj_met)

    # Year, Month, and Day as YYMMDD

    @property
    def ymd(self):
        """Returns the MET value as a string in the form of YYMMDD in UTC"""
        return self.datetime.strftime("%y%m%d")

    @classmethod
    def from_ymd(cls, ymd):
        """Create a Met object with its value set to the equivalent of YYMMDD string"""
        dt = datetime.datetime.strptime(ymd, '%y%m%d')
        return cls.from_datetime(dt)

    # Year, Month, Day, and Hour as YYMMDD_HH

    @property
    def ymd_h(self):
        """Returns the MET value as a string in the form of YYMMDD_HHz in UTC"""
        return self.datetime.strftime("%y%m%d_%Hz")

    @classmethod
    def from_ymd_h(cls, ymd):
        """Create a Met object with its value set to the equivalent of YYMMDD_HHz string"""
        dt = datetime.datetime.strptime(ymd, '%y%m%d_%Hz')
        return cls.from_datetime(dt)

    # Current time

    @classmethod
    def now(cls):
        """Create a Met object with its value set to the current time"""
        m = cls(0)
        m.__time = Time.now()
        return m

    # String functions

    def iso(self):
        """Returns the MET value as a string in the form of yyyy-mm-ddTHH:MM:SS in UTC"""
        return self.datetime.strftime("%Y-%m-%dT%H:%M:%S")

    def __repr__(self):
        """Returns a string representation of the Met object"""
        return "<Met seconds = {:.6f}>".format(self.met)

    # Math functions

    def add(self, x):
        """Returns an Met object with its value set to this object's value with x seconds added to it"""
        if not(isinstance(x, int) or isinstance(x, float)):
            raise ValueError("Can only add int or float to Met")
        return Met(self.met + x)

    def sub(self, x):
        """Returns an Met object with its value set to this object's value with x seconds subtracted from it"""
        if isinstance(x, Met):
            return self.met - x.met
        elif isinstance(x, int) or isinstance(x, float):
            return Met(self.met - x)
        raise ValueError("Can only subtract int, float or Met from Met")

    # Overriding built-in operators

    def __add__(self, other):
        return self.add(other)

    def __sub__(self, other):
        return self.sub(other)

    def __lt__(self, other):
        if isinstance(other, Met):
            return self.met < other.met
        else:
            raise TypeError("'<' not supported between instances of 'Met' and '{}'".format(type(other)))

    def __le__(self, other):
        if isinstance(other, Met):
            return self.met <= other.met
        else:
            raise TypeError("'<=' not supported between instances of 'Met' and '{}'".format(type(other)))

    def __gt__(self, other):
        if isinstance(other, Met):
            return self.met > other.met
        else:
            raise TypeError("'>' not supported between instances of 'Met' and '{}'".format(type(other)))

    def __ge__(self, other):
        if isinstance(other, Met):
            return self.met >= other.met
        else:
            raise TypeError("'>=' not supported between instances of 'Met' and '{}'".format(type(other)))

    def __eq__(self, other):
        if isinstance(other, Met):
            return self.met == other.met
        else:
            raise TypeError("'==' not supported between instances of 'Met' and '{}'".format(type(other)))

    def __ne__(self, other):
        if isinstance(other, Met):
            return self.met != other.met
        else:
            raise TypeError("'!=' not supported between instances of 'Met' and '{}'".format(type(other)))


# Some time related functions

def inclusive_date_range(start, stop, step=datetime.timedelta(days=1)):
    """Creates a list of Met from start to stop
    :param start:
    :param stop:
    :param step:
    :return:
    """
    d = start
    result = []

    if start <= stop:
        earliest, latest = start, stop
    else:
        earliest, latest = stop, start

    while earliest <= d <= latest:
        result.append(d)
        d += step

    return result


def dates_range_from(num_days, dt=datetime.datetime.utcnow().date()):
    """Creates a list of dates within the given range

    :param num_days: Number of days to include in the list
    :param dt: The last date to be included in the list
    :return: List of date values representing hours.
    """
    d = dt - datetime.timedelta(days=num_days - 1)
    return inclusive_date_range(d, dt)


def hours_range_from(num_hours, dt=datetime.datetime.utcnow()):
    """Creates a list of datetimes within the given range

    :param num_hours: Number of hours to include in the list
    :param dt: The last hour to be included in the list (datetime will be truncated to hour value)
    :return: List of datetime values representing hours.
    """
    d = datetime.datetime(dt.year, dt.month, dt.day, dt.hour, 0, 0)
    d -= datetime.timedelta(hours=num_hours - 1)

    return inclusive_date_range(d, dt, datetime.timedelta(hours=1))


def dates_from_hours(hours):
    """Converts a list of hours to a list of days spanned
    :param hours:
    :return: list of dates
    """
    return inclusive_date_range(hours[0].date(), hours[-1].date())
