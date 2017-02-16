# timing functions

current_milli_time = lambda: int(round(time.time() * 1000))

def milli_to_hr(millis):
    return (1.0 * millis) / (3600000.0)

def t_diff(start, finish, ts, steps):
    elapsed = milli_to_hr(finish-start)
    ns = (ts * steps * 1.0) / (1000000 * 1.0)
    return elapsed, elapsed / ns
