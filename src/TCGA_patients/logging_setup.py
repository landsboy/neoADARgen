import logging
from rich.logging import RichHandler

def setup_logging(verbose: bool = False, log_file: str = "logs/TCGA_patients.log") -> None:
    level = logging.DEBUG if verbose else logging.INFO

    handlers = [RichHandler(rich_tracebacks=True)]
    file_handler = logging.FileHandler(log_file, mode="w", encoding="utf-8")
    file_handler.setFormatter(logging.Formatter("%(asctime)s | %(levelname)-8s | %(name)s | %(message)s"))
    handlers.append(file_handler)

    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
        datefmt="%y-%m-%d %H:%M:%S",
        handlers=handlers,
    )
