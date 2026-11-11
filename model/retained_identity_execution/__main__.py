"""Run every enumerated original case with independent current validation."""
from .live import authority
from .runner import Execution, case_failure


def main():
    value = authority()
    nodes = value['baseline_collection_order']
    execution = Execution(nodes)
    results = execution.run()
    failures = {node: case_failure(phases) for node, phases in results.items() if case_failure(phases)}
    print(f'{len(nodes)} actual retained cases at {value["revision"]}; {len(failures)} failures')
    print(f'Historical scope; current dependency validation separate. Journal: {execution.log_directory}')
    if failures:
        for node, diagnostic in failures.items():
            print(node, diagnostic)
        raise SystemExit(1)


if __name__ == '__main__':
    main()
