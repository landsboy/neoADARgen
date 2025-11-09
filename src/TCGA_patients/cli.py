from pathlib import Path
import typer
from rich import print
from .config import AppConfig, load_config
from .logging_setup import setup_logging
from .io_utils import setup_environment
from .pipeline import main

app = typer.Typer(add_completion=False, no_args_is_help=True)

@app.command()
def run(
    config: Path = typer.Option(
        None,
        "--config", "-c",
        help="Path to YAML configuration file."
    ),
    project_dir: Path = typer.Option(
        None, "--project_dir", "-i", help="Override: folder with patient data."
    ),
    results_dir: Path = typer.Option(
        None, "--results_dir", "-o", help="Override: output folder for results."
    ),
    sup_dir: Path = typer.Option(
        None, "--sup_dir", "-s", help="Override: supplementary data directory."
    ),
    netmhc_path: Path = typer.Option(
        None, "--netmhc_path", "-p", help="Override: NetMHCpan installation path."
    ),
    hg38_fa: Path = typer.Option(
        None, "--hg38_fa", "-f", help="Override: hg38 FASTA reference file."
    ),
    num_editing: list[int] = typer.Option(
        None, "--num_editing", "-n", help="Override: edit modes (0/1/2)."
    ),
    verbose: bool = typer.Option(
        None, "--verbose", "-v", help="Override: verbose mode."
    ),
):
    """
    Run the pipeline using YAML config or optional CLI overrides.
    """
    # 1. Load from YAML or initialize defaults
    if config and config.exists():
        print(f"[yellow]Loading configuration from:[/yellow] {config}")
        cfg = load_config(config)
    else:
        print(f"[yellow]Using CLI arguments only.[/yellow]")
        cfg = AppConfig()

    # 2. Apply CLI overrides (only if provided)
    if project_dir is not None:
        cfg.paths.project_dir = project_dir
    if results_dir is not None:
        cfg.paths.results_dir = results_dir
    if sup_dir is not None:
        cfg.paths.sup_dir = sup_dir
    if netmhc_path is not None:
        cfg.paths.netmhc_path = netmhc_path
    if hg38_fa is not None:
        cfg.paths.hg38_fa = hg38_fa
    if num_editing is not None:
        cfg.runtime.edit_modes = num_editing
    if verbose is not None:
        cfg.runtime.verbose = verbose

    # Check that all required paths are set
    required_paths = {
        "project_dir": cfg.paths.project_dir,
        "results_dir": cfg.paths.results_dir,
        "sup_dir": cfg.paths.sup_dir,
        "netmhc_path": cfg.paths.netmhc_path,
        "hg38_fa": cfg.paths.hg38_fa,
    }

    missing = [name for name, path in required_paths.items() if path is None]
    if missing:
        print(f"[red]Missing required arguments:[/red] {', '.join(missing)}")
        print("[yellow]You must provide them either in config.yml or via CLI.[/yellow]")
        raise typer.Exit(code=1)
    
    # validate all important paths
    validate_paths(cfg)

    # 3. Logging & environment
    setup_logging(verbose=cfg.runtime.verbose, log_file=cfg.runtime.log_file)
    setup_environment(cfg.paths.netmhc_path)

    print(f"[bold cyan]Starting pipeline with config: {config}[/bold cyan]")
    main(cfg)
    print("[bold green]Pipeline completed successfully![/bold green]")

def validate_paths(cfg):
    errors = []

    if not cfg.paths.project_dir.exists():
        errors.append(f"[red] project_dir not found:[/red] {cfg.paths.project_dir}")

    if not cfg.paths.sup_dir.exists():
        errors.append(f"[red] sup_dir not found:[/red] {cfg.paths.sup_dir}")

    if not cfg.paths.netmhc_path.exists():
        errors.append(f"[red] NetMHCpan path not found:[/red] {cfg.paths.netmhc_path}")

    if not cfg.paths.hg38_fa.exists():
        errors.append(f"[red] hg38 FASTA file not found:[/red] {cfg.paths.hg38_fa}")

    # Create results_dir if missing
    if not cfg.paths.results_dir.exists():
        print(f"[yellow]⚠ results_dir not found, creating it at:[/yellow] {cfg.paths.results_dir}")
        cfg.paths.results_dir.mkdir(parents=True, exist_ok=True)

    if errors:
        for e in errors:
            print(e)
        raise typer.Exit(code=1)   # Stop the script if critical paths are missing

if __name__ == "__main__":
    app()
