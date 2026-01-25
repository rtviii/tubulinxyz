import asyncio
import os
from typing import Optional
import typer
from rich.console import Console
from neo4j import GraphDatabase, Driver
from typing_extensions import Annotated
from lib.etl.collector import TubulinETLCollector
from lib.etl.assets import GlobalOps, TubulinStructureAssets
from lib.etl.constants import NEO4J_URI, NEO4J_USER, NEO4J_CURRENTDB, NEO4J_PASSWORD
from neo4j_tubxz.db_lib_builder import Neo4jAdapter
from neo4j_tubxz.db_lib_reader import Neo4jReader
from rich.table import Table

from neo4j_tubxz.models import StructureFilters  # <--- ADD THIS IMPORT



app = typer.Typer(
    help="CLI for collecting and uploading Tubulin structure data.",
    add_completion=False
)
console = Console()  


def get_db_driver() -> Driver:
    """Initializes and returns a Neo4j Driver instance."""
    if not all([NEO4J_URI, NEO4J_USER, NEO4J_PASSWORD, NEO4J_CURRENTDB]):
        raise EnvironmentError("Missing Neo4j environment variables!")
    try:
        return GraphDatabase.driver(
            NEO4J_URI, 
            auth=(NEO4J_USER, NEO4J_PASSWORD), 
            database=NEO4J_CURRENTDB
        )
    except Exception as e:
        console.print(f"[bold red]Failed to create Neo4j driver:[/bold red] {e}")
        raise typer.Exit(code=1)

@app.command(name="init-db")
def init_db_command():
    """
    Initializes the Neo4j database: sets constraints and builds phylogenies.
    """
    try:
        console.print("Initializing Neo4j database (setting constraints, building phylogeny)...", style="cyan")
        adapter = get_db_adapter()
        adapter.initialize_new_instance() # This runs constraints AND phylogenies
        console.print("[bold green]Database initialized successfully.[/bold green]")
    except Exception as e:
        console.print(f"[bold red]Error initializing database:[/bold red] {e}")
        raise typer.Exit(code=1)

@app.command(name="collect-one")
def collect_one(
    rcsb_id: Annotated[str, typer.Argument(help="The 4-character PDB ID to collect.")],
    overwrite: Annotated[bool, typer.Option("--overwrite", "-o", help="Overwrite existing profile on disk.")] = False
):
    """
    Fetches a single structure from RCSB and saves its profile to disk.
    """
    rcsb_id = rcsb_id.upper()
    console.print(f"Attempting to collect [bold magenta]{rcsb_id}[/bold magenta]...", style="cyan")

    try:
        collector = TubulinETLCollector(rcsb_id)
        # Check if it exists *before* running the async collector
        if not overwrite and os.path.exists(collector.assets.paths.profile):
            console.print(f"[yellow]Skipped:[/yellow] Profile for {rcsb_id} already exists. Use -o to overwrite.")
            return

        profile = asyncio.run(collector.generate_profile(overwrite=overwrite))
        # console.print(f"[bold green]Successfully collected and saved profile:[/bold green] {profile.rcsb_id}")
        console.print(f"Profile saved at: {collector.assets.paths.profile}")
    
    except Exception as e:
        console.print(f"[bold red]Error collecting {rcsb_id}:[/bold red] {e}")
        import traceback
        traceback.print_exc()
        raise typer.Exit(code=1)

@app.command(name="collect-missing")
def collect_missing():
    """
    Finds all tubulin structures on RCSB that are missing locally,
    then collects and saves them all.
    """
    console.print("Finding missing profiles (RCSB vs. Local Disk)...", style="cyan")
    
    try:
        missing_ids = GlobalOps.missing_profiles()
    except Exception as e:
        console.print(f"[bold red]Error fetching list of missing profiles:[/bold red] {e}")
        raise typer.Exit(code=1)

    if not missing_ids:
        console.print("✨ [bold green]All local profiles are up-to-date with RCSB![/bold green]")
        return

    console.print(f"Found [bold yellow]{len(missing_ids)}[/bold yellow] new structures to collect.")
    
    success_count = 0
    fail_count = 0
    failed_ids = []

    for rcsb_id in missing_ids:
        console.print(f"--- Collecting [bold magenta]{rcsb_id}[/bold magenta]... ---")
        try:
            collector = TubulinETLCollector(rcsb_id)
            asyncio.run(collector.generate_profile(overwrite=False))
            console.print(f"[green]Successfully collected {rcsb_id}[/green]")
            success_count += 1
        except Exception as e:
            console.print(f"[bold red]!!! Failed to collect {rcsb_id}:[/bold red] {e}")
            fail_count += 1
            failed_ids.append(rcsb_id)
    
    console.print("\n--- [bold]Collection Complete[/bold] ---")
    console.print(f"[green]Successful:[/green] {success_count}")
    console.print(f"[red]Failed:[/red] {fail_count}")
    if failed_ids:
        console.print(f"[red]Failed IDs:[/red] {', '.join(failed_ids)}")

def get_all_structs_in_db(driver: Driver) -> list[str]:
    """
    Fetches a list of all rcsb_id strings from Structure nodes in Neo4j.
    """
    console.print("Querying database for existing structure IDs...", style="cyan")
    query = "MATCH (s:Structure) RETURN s.rcsb_id AS rcsb_id"
    
    try:
        with driver.session() as session:
            results = session.run(query)
            struct_ids = [record["rcsb_id"] for record in results]
            console.print(f"Found [bold yellow]{len(struct_ids)}[/bold yellow] structures in database.")
            return struct_ids
    except Exception as e:
        console.print(f"[bold red]Error querying database:[/bold red] {e}")
        raise typer.Exit(code=1)


from concurrent.futures import ThreadPoolExecutor, as_completed

@app.command(name="upload-missing")
def upload_missing(
    workers: Annotated[int, typer.Option("--workers", "-w", help="Number of parallel workers")] = 8
):
    """
    Finds all structures on disk that are missing from the Neo4j database,
    then uploads them all in parallel.
    """
    console.print("Finding local profiles missing from the database...", style="cyan")
    db_driver = None
    
    try:
        local_profiles = set(GlobalOps.list_profiles())
        if not local_profiles:
            console.print("[yellow]No local profiles found in TUBETL_DATA.[/yellow]")
            return

        db_driver = get_db_driver()
        db_profiles = set(get_all_structs_in_db(db_driver))
        
        missing_from_db = sorted(list(local_profiles - db_profiles))

        if not missing_from_db:
            console.print("[bold green]Database is fully synced with local profiles![/bold green]")
            return

        console.print(f"Found [bold yellow]{len(missing_from_db)}[/bold yellow] structures to upload with {workers} workers.")
        
        success_count = 0
        fail_count = 0
        failed_ids = []

        def upload_one_structure(rcsb_id: str) -> tuple[str, bool, str]:
            """Returns (rcsb_id, success, error_msg)"""
            try:
                adapter = get_db_adapter()  # Each thread gets its own adapter
                adapter.add_total_structure(rcsb_id, disable_exists_check=True)
                return (rcsb_id, True, "")
            except Exception as e:
                return (rcsb_id, False, str(e))

        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(upload_one_structure, rid): rid for rid in missing_from_db}
            
            for future in as_completed(futures):
                rcsb_id, success, error_msg = future.result()
                if success:
                    console.print(f"[green]Uploaded {rcsb_id}[/green]")
                    success_count += 1
                else:
                    console.print(f"[red]Failed {rcsb_id}:[/red]", markup=True)
                    console.print(error_msg, markup=False)
                    fail_count += 1
                    failed_ids.append(rcsb_id)

        console.print("\n--- [bold]Upload Complete[/bold] ---")
        console.print(f"[green]Successful:[/green] {success_count}")
        console.print(f"[red]Failed:[/red] {fail_count}")
        if failed_ids:
            console.print(f"[red]Failed IDs:[/red] {', '.join(failed_ids)}")

    except Exception as e:
        console.print("[bold red]An error occurred during the upload-missing operation:[/bold red]")
        console.print(str(e), markup=False)
        import traceback
        traceback.print_exc()
        raise typer.Exit(code=1)
    finally:
        if db_driver:
            db_driver.close()
            console.print("Database connection closed.", style="dim")

@app.command(name="db-stats")
def db_stats():
    """
    Display basic statistics about the Neo4j database content.
    """
    reader = Neo4jReader()
    try:
        console.print("Fetching database statistics...", style="cyan")
        
        # Node Labels
        labels = reader.node_types()
        console.print(f"[bold]Node Labels:[/bold] {', '.join(labels)}")
        
        # ID Count
        ids = reader.all_ids()
        console.print(f"[bold]Total Structures:[/bold] {len(ids)}")
        if ids:
            console.print(f"[dim]Sample IDs: {', '.join(ids[:5])}...[/dim]")
            
        # Taxa Counts
        src_taxa = reader.get_taxa("source")
        host_taxa = reader.get_taxa("host")
        console.print(f"[bold]Unique Source Taxa:[/bold] {len(src_taxa)}")
        console.print(f"[bold]Unique Host Taxa:[/bold] {len(host_taxa)}")

    except Exception as e:
        console.print(f"[bold red]Error fetching stats:[/bold red] {e}")

@app.command(name="inspect-struct")
def inspect_structure(
    rcsb_id: Annotated[Optional[str], typer.Argument(help="RCSB ID to inspect. If omitted, picks random.")] = None
):
    """
    Fetch a structure (random or specific) and display its deep data (ligands, proteins).
    """
    reader = Neo4jReader()
    try:
        if rcsb_id:
            console.print(f"Fetching [magenta]{rcsb_id}[/magenta]...", style="cyan")
            data = reader.get_structure_by_id(rcsb_id)
            if not data:
                console.print(f"[red]Structure {rcsb_id} not found in DB.[/red]")
                return
        else:
            console.print("Fetching [magenta]random structure[/magenta]...", style="cyan")
            data = reader.random_structure()
        
        # Display
        console.print(f"\n[bold underline]Structure: {data['rcsb_id']}[/bold underline]")
        console.print(f"Title: {data.get('citation_title', 'N/A')}")
        console.print(f"Method: {data.get('expMethod', 'N/A')} | Res: {data.get('resolution', 'N/A')}A")
        
        proteins = data.get('proteins', [])
        console.print(f"\n[bold]Proteins ({len(proteins)}):[/bold]")
        for p in proteins:
            fam = p.get('family', 'Unknown')
            console.print(f" - Chain [yellow]{p['auth_asym_id']}[/yellow]: {p.get('rcsb_pdbx_description', 'N/A')} (Family: {fam})")

        ligands = data.get('ligands', [])
        console.print(f"\n[bold]Ligands ({len(ligands)}):[/bold]")
        for l in ligands:
            console.print(f" - [green]{l['chemicalId']}[/green]: {l.get('chemicalName', 'N/A')}")

    except Exception as e:
        console.print(f"[bold red]Error:[/bold red] {e}")

@app.command(name="list-structs")
def list_structures(
    limit: int = 20,
    search: Optional[str] = None,
    min_year: Optional[int] = None,
    source_taxid: Optional[int] = None
):
    """
    Search structures with filters.
    """
    reader = Neo4jReader()
    
    # Build params
    filters = StructureFilters(
        limit=limit,
        search=search,
        year_min=min_year if min_year else None,
        sourceTaxa=[source_taxid] if source_taxid else None
    )
    
    try:
        structs, total, next_cursor = reader.list_structs_filtered(filters)
        console.print(f"Found [bold]{total}[/bold] structures matching criteria.")
        
        # --- FIX IS HERE: Use Table directly ---
        table = Table("ID", "Year", "Resolution", "Title")
        
        for s in structs:
            table.add_row(
                s['rcsb_id'], 
                str(s.get('citation_year', '')), 
                str(s.get('resolution', '')), 
                (s.get('citation_title', 'N/A') or "N/A")[:50]+"..."
            )
        console.print(table)
        
        if next_cursor:
            console.print(f"[dim]Next cursor: {next_cursor}[/dim]")
            
    except Exception as e:
        console.print(f"[bold red]Error:[/bold red] {e}")

@app.command(name="collect-all")
def collect_all(
    overwrite: Annotated[bool, typer.Option("--overwrite", "-o", help="Overwrite existing profiles on disk.")] = False
):
    """
    Collects ALL tubulin structures from RCSB, optionally overwriting existing profiles.
    """
    console.print("Fetching complete list of tubulin structures from RCSB...", style="cyan")
    
    try:
        all_rcsb_ids = GlobalOps.current_rcsb_structs()
    except Exception as e:
        console.print(f"[bold red]Error fetching RCSB structure list:[/bold red] {e}")
        raise typer.Exit(code=1)

    if not all_rcsb_ids:
        console.print("[yellow]No tubulin structures found on RCSB.[/yellow]")
        return

    console.print(f"Found [bold yellow]{len(all_rcsb_ids)}[/bold yellow] total structures on RCSB.")
    
    # Filter out existing if not overwriting
    to_process = all_rcsb_ids
    if not overwrite:
        existing = set(GlobalOps.list_profiles())
        to_process = [rid for rid in all_rcsb_ids if rid not in existing]
        console.print(f"Skipping {len(existing)} existing profiles. Processing [bold yellow]{len(to_process)}[/bold yellow] structures.")
    else:
        console.print(f"Overwrite enabled. Processing all [bold yellow]{len(to_process)}[/bold yellow] structures.")

    if not to_process:
        console.print("✨ [bold green]All profiles already exist![/bold green]")
        return

    success_count = 0
    fail_count = 0
    failed_ids = []

    for i, rcsb_id in enumerate(to_process, 1):
        console.print(f"\n[{i}/{len(to_process)}] Collecting [bold magenta]{rcsb_id}[/bold magenta]...")
        try:
            collector = TubulinETLCollector(rcsb_id)
            asyncio.run(collector.generate_profile(overwrite=overwrite))
            console.print(f"[green]✓ {rcsb_id}[/green]")
            success_count += 1
        except Exception as e:
            console.print(f"[bold red]✗ Failed {rcsb_id}:[/bold red] {e}")
            fail_count += 1
            failed_ids.append(rcsb_id)
    
    console.print("\n--- [bold]Collection Complete[/bold] ---")
    console.print(f"[green]Successful:[/green] {success_count}")
    console.print(f"[red]Failed:[/red] {fail_count}")
    if failed_ids:
        console.print(f"[red]Failed IDs:[/red] {', '.join(failed_ids)}")



# In cli.py

from lib.etl.constants import NEO4J_CURRENTDB, NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER
from neo4j_tubxz.db_lib_builder import Neo4jAdapter
from lib.etl.assets import GlobalOps, TubulinStructureAssets
import os
def get_db_adapter() -> Neo4jAdapter:
    """Factory for getting the Neo4j adapter."""
    return Neo4jAdapter(NEO4J_URI, NEO4J_USER, NEO4J_CURRENTDB, NEO4J_PASSWORD)


@app.command(name="init-db")
def init_db(
    db_name: Annotated[str, typer.Option("--name", "-n", help="Name for new database (default: use NEO4J_CURRENTDB from constants)")] = "",
    skip_phylogeny: Annotated[bool, typer.Option("--skip-phylogeny", help="Skip phylogeny seeding (faster, but no taxonomy filtering)")] = False,
    create: Annotated[bool, typer.Option("--create", "-c", help="Create the database if it doesn't exist")] = False,
):
    """
    Initialize a Neo4j database with constraints and phylogeny tree.
    
    Use --create to create a new database first (requires Enterprise/Aura).
    """
    from neo4j import GraphDatabase
    
    target_db = db_name if db_name else NEO4J_CURRENTDB
    
    console.print(f"Target database: [bold]{target_db}[/bold] at [bold]{NEO4J_URI}[/bold]", style="cyan")
    
    try:
        # If --create flag, create the database first via system db
        if create:
            console.print(f"Creating database [bold]{target_db}[/bold]...")
            
            # Connect to system database to create new db
            system_driver = GraphDatabase.driver(
                NEO4J_URI, 
                auth=(NEO4J_USER, NEO4J_PASSWORD), 
                database="system"
            )
            
            with system_driver.session() as session:
                # Check if database already exists
                result = session.run("SHOW DATABASES YIELD name")
                existing_dbs = [r["name"] for r in result]
                
                if target_db in existing_dbs:
                    console.print(f"[yellow]Database '{target_db}' already exists.[/yellow]")
                else:
                    session.run(f"CREATE DATABASE {target_db}")
                    console.print(f"[green]Created database '{target_db}'.[/green]")
                    
                    # Wait for database to come online
                    console.print("Waiting for database to start...")
                    import time
                    for _ in range(30):  # Wait up to 30 seconds
                        result = session.run(
                            "SHOW DATABASE $name YIELD currentStatus", 
                            {"name": target_db}
                        ).single()
                        if result and result["currentStatus"] == "online":
                            break
                        time.sleep(1)
                    else:
                        console.print("[yellow]Warning: Database may still be starting up.[/yellow]")
            
            system_driver.close()
        
        # Now connect to the target database and initialize
        adapter = Neo4jAdapter(NEO4J_URI, NEO4J_USER, target_db, NEO4J_PASSWORD)
        
        console.print("Creating constraints and indexes...")
        adapter.init_constraints()
        console.print("[green]Constraints created.[/green]")
        
        if not skip_phylogeny:
            console.print("Seeding phylogeny tree (this may take a minute)...")
            adapter.init_phylogenies()
            console.print("[green]Phylogeny tree seeded.[/green]")
        else:
            console.print("[yellow]Skipped phylogeny seeding.[/yellow]")
        
        adapter.close()
        console.print(f"[bold green]Database '{target_db}' initialized successfully![/bold green]")
        
        # Remind user to update constants if they used a custom name
        if db_name and db_name != NEO4J_CURRENTDB:
            console.print(f"\n[yellow]Note: Update NEO4J_CURRENTDB in lib/etl/constants.py to '{db_name}' to use this database by default.[/yellow]")
        
    except Exception as e:
        console.print(f"[bold red]Error:[/bold red] {e}")
        import traceback
        traceback.print_exc()
        raise typer.Exit(code=1)


@app.command(name="upload-one")
def upload_one(
    rcsb_id: Annotated[str, typer.Argument(help="The 4-character PDB ID to upload from disk.")],
    overwrite: Annotated[bool, typer.Option("--overwrite", "-o", help="Overwrite if structure already exists")] = False,
):
    """
    Upload a single structure profile to Neo4j.
    """
    rcsb_id = rcsb_id.upper()
    console.print(f"Uploading [bold magenta]{rcsb_id}[/bold magenta] to Neo4j...", style="cyan")

    # Verify profile exists
    try:
        assets = TubulinStructureAssets(rcsb_id)
        if not os.path.exists(assets.paths.profile):
            raise FileNotFoundError(f"Profile not found at {assets.paths.profile}")
        console.print(f"Found profile: {assets.paths.profile}")
    except FileNotFoundError as e:
        console.print(f"[bold red]Error:[/bold red] {e}")
        raise typer.Exit(code=1)

    # Upload
    try:
        adapter = get_db_adapter()
        adapter.add_total_structure(rcsb_id, disable_exists_check=overwrite)
        adapter.close()
        console.print(f"[bold green]Successfully uploaded {rcsb_id}.[/bold green]")
    except Exception as e:
        console.print(f"[bold red]Error uploading {rcsb_id}:[/bold red] {e}")
        import traceback
        traceback.print_exc()
        raise typer.Exit(code=1)


@app.command(name="upload-all")
def upload_all(
    overwrite: Annotated[bool, typer.Option("--overwrite", "-o", help="Overwrite existing structures")] = False,
    limit: Annotated[int, typer.Option("--limit", "-n", help="Limit number of structures to upload")] = 0,
):
    """
    Upload all structure profiles from disk to Neo4j.
    """
    profiles = GlobalOps.list_profiles()
    
    if limit > 0:
        profiles = profiles[:limit]
    
    console.print(f"Found [bold]{len(profiles)}[/bold] profiles to upload.", style="cyan")
    
    if not profiles:
        console.print("[yellow]No profiles found on disk.[/yellow]")
        raise typer.Exit(code=0)

    adapter = get_db_adapter()
    success_count = 0
    error_count = 0

    for i, rcsb_id in enumerate(sorted(profiles), 1):
        console.print(f"[{i}/{len(profiles)}] {rcsb_id}...", end=" ")
        try:
            adapter.add_total_structure(rcsb_id, disable_exists_check=overwrite)
            console.print("[green]OK[/green]")
            success_count += 1
        except Exception as e:
            console.print(f"[red]FAILED: {e}[/red]")
            error_count += 1

    adapter.close()
    console.print(f"\n[bold]Done:[/bold] {success_count} succeeded, {error_count} failed.")


@app.command(name="delete-one")
def delete_one(
    rcsb_id: Annotated[str, typer.Argument(help="The 4-character PDB ID to delete.")],
    force: Annotated[bool, typer.Option("--force", "-f", help="Skip confirmation")] = False,
):
    """
    Delete a structure from Neo4j.
    """
    rcsb_id = rcsb_id.upper()
    
    if not force:
        confirm = typer.confirm(f"Are you sure you want to delete {rcsb_id}?")
        if not confirm:
            raise typer.Abort()

    try:
        adapter = get_db_adapter()
        if not adapter.check_structure_exists(rcsb_id):
            console.print(f"[yellow]Structure {rcsb_id} not found in database.[/yellow]")
            raise typer.Exit(code=0)
        
        adapter.delete_structure(rcsb_id)
        adapter.close()
        console.print(f"[bold green]Deleted {rcsb_id}.[/bold green]")
    except Exception as e:
        console.print(f"[bold red]Error:[/bold red] {e}")
        raise typer.Exit(code=1)


@app.command(name="db-status")
def db_status():
    """
    Show database connection status and basic stats.
    """
    console.print(f"Connecting to [bold]{NEO4J_URI}[/bold] / [bold]{NEO4J_CURRENTDB}[/bold]...", style="cyan")
    
    try:
        adapter = get_db_adapter()
        
        # Test connection and get counts
        with adapter.driver.session() as session:
            result = session.run("""
                MATCH (s:Structure) WITH count(s) AS structures
                MATCH (e:Entity) WITH structures, count(e) AS entities
                MATCH (i:Instance) WITH structures, entities, count(i) AS instances
                MATCH (c:Chemical) WITH structures, entities, instances, count(c) AS chemicals
                MATCH (p:PhylogenyNode) WITH structures, entities, instances, chemicals, count(p) AS taxa
                MATCH (v:Variant) WITH structures, entities, instances, chemicals, taxa, count(v) AS variants
                RETURN structures, entities, instances, chemicals, taxa, variants
            """).single()
            
            console.print("\n[bold green]Connected successfully![/bold green]\n")
            console.print(f"  Structures:  {result['structures']}")
            console.print(f"  Entities:    {result['entities']}")
            console.print(f"  Instances:   {result['instances']}")
            console.print(f"  Chemicals:   {result['chemicals']}")
            console.print(f"  Taxa:        {result['taxa']}")
            console.print(f"  Variants:    {result['variants']}")
        
        adapter.close()
        
    except Exception as e:
        console.print(f"[bold red]Connection failed:[/bold red] {e}")
        raise typer.Exit(code=1)


if __name__ == "__main__":
    # Ensure environment variables are set
    if not os.environ.get("TUBETL_DATA") or not os.environ.get("NEO4J_URI"):
        console.print("[bold red]Error: Environment variables (e.g., TUBETL_DATA, NEO4J_URI) are not set.[/bold red]")
        console.print("Please source your environment file before running the CLI.")
        exit(1)
    app()
