#!/bin/bash

# FAST TCP CLI Test Script
# This script provides convenient testing commands for the CLI

set -e  # Exit on any error

# Activate virtual environment if it exists
if [ -d "virtual" ]; then
    echo "Activating virtual environment..."
    source virtual/bin/activate
fi

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}FAST TCP CLI Test Suite${NC}"
echo "=========================="

# Function to print colored output
print_status() {
    if [ $1 -eq 0 ]; then
        echo -e "${GREEN}✓ $2${NC}"
    else
        echo -e "${RED}✗ $2${NC}"
    fi
}

# Check if Python environment is set up
echo "Checking Python environment..."
if ! python -c "import fast_tcp" 2>/dev/null; then
    echo -e "${YELLOW}Warning: fast_tcp module not found. Installing in development mode...${NC}"
    pip install -e .
fi

# Handle different test modes
if [ "$1" = "--quick" ]; then
    echo -e "\n${YELLOW}Running quick tests only...${NC}"
    if python quick_test.py; then
        print_status 0 "Quick tests completed successfully"
    else
        print_status 1 "Quick tests failed"
        exit 1
    fi
elif [ "$1" = "--all" ] || [ "$1" = "--comprehensive" ]; then
    echo -e "\n${YELLOW}Running comprehensive algorithm tests...${NC}"
    if python test_cli.py --all-algorithms; then
        print_status 0 "Comprehensive tests completed successfully"
    else
        print_status 1 "Comprehensive tests failed"
        exit 1
    fi
else
    # Run standard tests by default
    echo -e "\n${YELLOW}Running standard CLI tests...${NC}"
    if python test_cli.py; then
        print_status 0 "Standard tests completed successfully"
    else
        print_status 1 "Standard tests failed"
        exit 1
    fi
fi



# Manual testing section
if [ "$1" = "--manual" ]; then
    echo -e "\n${YELLOW}Manual Testing Mode${NC}"
    echo "You can now test the CLI manually. Here are some examples:"
    echo ""
    echo "1. Test help:"
    echo "   python -m fast_tcp.cli --help"
    echo ""
    echo "2. Test with sample data:"
    echo "   python -m fast_tcp.cli --test-dir input/chart_v0 --algo FAST-pw --entity function --repetitions 1"
    echo ""
    echo "3. Test with custom output:"
    echo "   python -m fast_tcp.cli --test-dir input/chart_v0 --algo STR --entity bbox --output-dir ./test_output"
    echo ""
    echo "4. Available test modes:"
    echo "   ./test_cli.sh --quick         # Fast tests (~10 seconds)"
    echo "   ./test_cli.sh                 # Standard tests (~1 minute)" 
    echo "   ./test_cli.sh --all           # Comprehensive tests (~5-10 minutes)"
    echo ""
    echo "Press Ctrl+C to exit manual mode."
    echo ""
    
    # Keep shell open for manual testing
    bash
fi

echo -e "\n${GREEN}All tests completed successfully!${NC}"