#!/usr/bin/env python3

import sys
from math import sqrt

def f(x,y,t):
    r = sqrt(x**2 + y**2)
    return 2**(8/3)/(t**(4/3)*(1+2*(t**2+r**2)+(t**2-r**2)**2)**(4/3))


def main() -> int:
    for l in sys.stdin:
        [x,y,t] = map(lambda x: float(x), l.split(" "))
        print(f(x,y,t))
    return 0


if __name__ == '__main__':
    sys.exit(main())
