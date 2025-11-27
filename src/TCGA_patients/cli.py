from pathlib import Path
import typer
from rich import print
from .config import AppConfig, load_config
from .logging_setup import setup_logging
from .io_utils import setup_environment
from .pipeline import main
from types import SimpleNamespace

app = typer.Typer(add_completion=False, no_args_is_help=True)

def validate_paths(cfg):
    errors = []

    if cfg.mode == "tcga":
        if not cfg.paths.project_dir.exists():
            errors.append(f"[red] project_dir not found:[/red] {cfg.paths.project_dir}")

        if not cfg.paths.sup_dir.exists():
            errors.append(f"[red] sup_dir not found:[/red] {cfg.paths.sup_dir}")

    if not cfg.paths.netmhc_path.exists():
        errors.append(f"[red] NetMHCpan path not found:[/red] {cfg.paths.netmhc_path}")

    if not cfg.paths.hg38_fa.exists():
        errors.append(f"[red] hg38 FASTA file not found:[/red] {cfg.paths.hg38_fa}")

    # results_dir: create if missing
    if not cfg.paths.results_dir.exists():
        print(f"[yellow]⚠ results_dir not found, creating it:[/yellow] {cfg.paths.results_dir}")
        cfg.paths.results_dir.mkdir(parents=True, exist_ok=True)

    # Fail on errors
    if errors:
        for e in errors:
            print(e)
        raise typer.Exit(code=1)


@app.command()
def run(
    config: Path = typer.Option(
        None, "--config", "-c", help="Path to YAML configuration file."
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
    mutation: str = typer.Option(
        None, "--mutation", "-m",
        help="Mutation(s) in format like chr1:g.12825232T>G"
    ),
    hla: list[str] = typer.Option(
    None, "--hla",
    help="Patient-specific HLA alleles.",
    ),
    
    gene_counts: Path = typer.Option(
        None, "--gene_counts",
        help="Optional gene-count file (TSV). If omitted, TPM step is skipped."
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
    def override(target, name, value):
        if value is not None:
            setattr(target, name, value)

    override(cfg.paths, "project_dir", project_dir)
    override(cfg.paths, "results_dir", results_dir)
    override(cfg.paths, "sup_dir", sup_dir)
    override(cfg.paths, "netmhc_path", netmhc_path)
    override(cfg.paths, "hg38_fa", hg38_fa)

    override(cfg.runtime, "edit_modes", num_editing)
    override(cfg.runtime, "verbose", verbose)

    # 3. Determine execution mode
    single_mode = mutation or hla or gene_counts
    cfg.mode = "single" if single_mode else "tcga"

    if cfg.mode == "single":
        print("[bold magenta]Running in SINGLE-MUTATION mode[/bold magenta]")

        # Required for single mode
        if mutation is None:
            print("[red]Missing required --mutation[/red]")
            raise typer.Exit(1)

        if not hla:
            print("[red]At least one --hla allele is required[/red]")
            raise typer.Exit(1)

        # Attach single-mode block
        cfg.single = SimpleNamespace(
            mutation=mutation,
            hla=hla,
            gene_counts=gene_counts
        )

        required_fields = ["results_dir","sup_dir", "netmhc_path", "hg38_fa"]

    else:
        print("[bold cyan]Running in TCGA mode[/bold cyan]")

        required_fields = ["project_dir", "results_dir", "sup_dir", "netmhc_path", "hg38_fa"]


    # 4. Check that all required paths are set
    missing = []
    for field in required_fields:
        if getattr(cfg.paths, field) is None:
            missing.append(field)

    if missing:
        print(f"[red]Missing required arguments:[/red] {', '.join(missing)}")
        raise typer.Exit(1)
    
    # validate all important paths
    validate_paths(cfg)

    # 5. Logging & environment
    setup_logging(verbose=cfg.runtime.verbose, log_file=cfg.runtime.log_file)
    setup_environment(cfg.paths.netmhc_path)

    print(f"[bold cyan]Starting pipeline:[/bold cyan]")
    main(cfg)
    print("[bold green]Pipeline completed successfully![/bold green]")

if __name__ == "__main__":
    app()
