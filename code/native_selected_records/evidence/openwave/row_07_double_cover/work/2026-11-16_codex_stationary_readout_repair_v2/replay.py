"""Complete read-only verification entry point; no creation options."""
import audit_saved
import verify_v6

if __name__=='__main__':
    audit_saved.main()
    verify_v6.main()
