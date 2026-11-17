#!/usr/bin/env python3
"""Read-only exact interface replay; provider files are never written."""
from verify import main
if __name__=='__main__':
    result=main()
    assert result['elapsed_seconds']<60,'headline exceeds the registered wall-time gate'
