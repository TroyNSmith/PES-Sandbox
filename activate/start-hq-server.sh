#!/usr/bin/env bash
set -e

SERVER_ROOT="$PWD/.server"
mkdir -p "$SERVER_ROOT"

if ! hq server status >/dev/null 2>&1; then
  echo "Starting HyperQueue server…"
  nohup hq server start --server-dir "$SERVER_ROOT" >/dev/null 2>&1 &
  sleep 1
fi

export HQ_SERVER_DIR="$SERVER_ROOT/hq-current"
