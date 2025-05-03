# Biomart MCP

[![smithery badge](https://smithery.ai/badge/@jzinno/biomart-mcp)](https://smithery.ai/server/@jzinno/biomart-mcp)
[![Verified on MseeP](https://mseep.ai/badge.svg)](https://mseep.ai/app/ad562201-95a7-4ceb-83d6-09eddc7fdac4)

### A MCP server to interface with Biomart

[Model Context Protocol](https://modelcontextprotocol.io/introduction) (MCP) is an open protocol that standardizes how applications provide context to LLMs developed by [Anthropic](https://www.anthropic.com/). Here we use the [MCP python-sdk](https://github.com/modelcontextprotocol/python-sdk) to create a MCP server that interfaces with Biomart via the [pybiomart](https://github.com/jrderuiter/pybiomart) package.

![Demo showing biomart-mcp in action](./assets/demo.png)

There is a short [demo video](assets/mcp-demo.mp4) showing the MCP server in action on Claude Desktop.

## Installation

### Installing via Smithery

To install Biomart MCP for Claude Desktop automatically via [Smithery](https://smithery.ai/server/@jzinno/biomart-mcp):

```bash
npx -y @smithery/cli install @jzinno/biomart-mcp --client claude
```

### Clone the repository

```bash
git clone https://github.com/jzinno/biomart-mcp.git
cd biomart-mcp
```

### Claude Desktop

```bash
uv run --with mcp[cli] mcp install --with pybiomart biomart-mcp.py
```

### Cursor

Via Cusror's agent mode, other models can take advantage of MCP servers as well, such as those form OpenAI or DeepSeek. Click the cursor setting cogwheel and naviagate to `MCP` and either add the MCP server to the global config or add it to the a project scope by adding `.cursor/mcp.json` to the project.

Example `.cursor/mcp.json`:

```json
{
    "mcpServers": {
        "Biomart": {
            "command": "uv",
            "args": [
                "run",
                "--with",
                "mcp[cli]",
                "--with",
                "pybiomart",
                "mcp",
                "run",
                "/your/path/to/biomart-mcp.py"
            ]
        }
    }
}
```

### Glama

<a href="https://glama.ai/mcp/servers/v5a3mlxviu">
  <img width="380" height="200" src="https://glama.ai/mcp/servers/v5a3mlxviu/badge" alt="Biomart MCP server" />
</a>

### Development

```bash
# Create a virtual environment
uv venv

# MacOS/Linux
source .venv/bin/activate

# Windows
.venv\Scripts\activate

uv sync #or uv add mcp[cli] pybiomart

# Run the server in dev mode
mcp dev biomart-mcp.py
```

## Features

Biomart-MCP provides several tools to interact with Biomart databases:

- **Mart and Dataset Discovery**: List available marts and datasets to explore the Biomart database structure
- **Attribute and Filter Exploration**: View common or all available attributes and filters for specific datasets
- **Data Retrieval**: Query Biomart with specific attributes and filters to get biological data
- **ID Translation**: Convert between different biological identifiers (e.g., gene symbols to Ensembl IDs)

## Contributing

Pull requests are welcome! Some small notes on development:

- We are only using `@mcp.tool()` here by design, this is to maximize compatibility with clients that support MCP as seen in the [docs](https://modelcontextprotocol.io/clients).
- We are using `@lru_cache` to cache results of functions that are computationally expensive or make external API calls.
- We need to be mindful to not blow up the context window of the model, for example you'll see `df.to_csv(index=False).replace("\r", "")` in many places. This csv style return is much more token efficient than something like `df.to_string()` where the majority of the tokens are whitespace. Also be mindful of the fact that pulling all genes from a chromosome or similar large request will also be too large for the context window.

## Potential Future Features

There of course many more features that could be added, some maybe beyond the scope of the name `biomart-mcp`. Here are some ideas:

- Add webscraping for resource sites with `bs4`, for example we got the Ensembl gene ID for NOTCH1 then maybe in some cases it would be usful to grap the collated `Comments and Description Text from UniProtKB` section from [it's page on UCSC](https://genome.ucsc.edu/cgi-bin/hgGene?db=hg38&hgg_chrom=chr9&hgg_gene=ENST00000651671.1&hgg_start=136494433&hgg_end=136546048&hgg_type=knownGene)
- $...$
