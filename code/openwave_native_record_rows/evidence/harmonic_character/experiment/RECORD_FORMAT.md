# Records

INITIALIZATION.json contains all16 new three-event native branch histories, initial/final unit positions, placement words, currents and decoded cyclic origins.

inputs/compatible_records.jsonl.gz is a byte-for-byte copy of the original complete source/control archive. Its rows are[mode, initial_unit_port, twelve_bit_word, twelve_currents, three_ternary_values, compatible]. The frozen hash and receipt preserve its earlier provenance and cost.

COUNTERS.bin.gz has131072 consecutive3-byte records, in archived-word order and origin0..7 order. The bytes are origin, final_counter, flags. Flags=2*compatible+int(phase_positive), where phase is the k4 boundary character. Incompatible words retain their counter and positive hold phase but are excluded from the conditional character average. Counter records are observer updates, not native endpoint calls.
