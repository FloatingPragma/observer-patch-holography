# Snapshot link context

CRITERIA_CAUTIONS.md and CRITERIA_SUBMISSION_CHECKLIST.md preserve the exact original bytes from openwave/. Their relative links retain that original repository-root context; they are not relative to this packet. The first packaging link check incorrectly treated those snapshot links as packet-local and stopped before writing a manifest. The corrected packaging check preserves the snapshots and resolves their links against openwave/. Scientific execution, receipts and gates were unaffected.
