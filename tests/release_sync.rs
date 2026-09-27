//! The hand-written numbers and the hand-written prose have to say what the code says.
//!
//! Everything checked here is a fact the compiler cannot see, which is exactly why it went stale:
//! `README.md` carried an install snippet and a docs.rs link pinned to `0.6.0` for two minor
//! versions, `mu_r` was documented with `t_forward_bond_vol`'s sentence (misspelled), a `*_now`
//! pricer promised a future date, a receiver swaption was labelled the payer side, and a coupon
//! schedule parameter was documented backwards from the kernel that reads it.  None of that breaks
//! a build.  All of it misleads a caller who is reading the docs *because* they cannot tell from
//! the signature, so all of it is cheap to assert and expensive to get wrong.
//!
//! The rules split into two groups.
//!
//! **Release sync** — the version string cargo reads from `package.version`, and the two places it is
//! written by hand (`README.md`'s install snippet and its docs.rs link), plus the CHANGELOG section
//! for the release and the empty `## [Unreleased]` block that the next one is written into.  See
//! the "Doc and version bookkeeping" section of `README.md` for the per-release sequence.
//!
//! **Doc-comment shape** — the copy/paste patterns that produced the drift above, checked
//! mechanically against every function in `src/` rather than against the handful that happened to be
//! noticed: a `*_now` function whose summary line talks about a future date, a receiver-side
//! function whose summary line says "payer" (and the reverse), and the misspelling "volality".
//!
//! The prose-length rule (`every_public_function_has_a_summary_that_says_something`) applies only to
//! the **public API** — the names `src/lib.rs` re-exports at the crate root and the modules declared
//! `pub mod`, which is the same set the `missing_docs` lint guards.  A `pub fn` inside a private
//! module (`mc`, `validation`, `schedules`' internals) is not `hull_white::` anything; its comment
//! is written for whoever edits that file next, and holding it to a character count would be noise
//! that drowns out the checks that are about a reader.
//!
//! Deliberately *not* checked here: that a public item has a doc comment at all.  That one is
//! enforced by the compiler (`[lints.rust]` in `Cargo.toml`) and does not need a test restating it.

use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::{Path, PathBuf};

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn read(path: &Path) -> String {
    fs::read_to_string(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()))
}

/// `package.version` from `Cargo.toml`, read as text rather than parsed as TOML: the manifest is a
/// dependency-free file to grep, and pulling in a TOML parser to police two lines in a README would
/// be a heavier thing to keep current than the drift being prevented.
fn crate_version() -> String {
    let manifest = read(&root().join("Cargo.toml"));
    let mut in_package = false;
    for line in manifest.lines() {
        let line = line.trim();
        if line.starts_with('[') {
            in_package = line == "[package]";
            continue;
        }
        if in_package {
            if let Some(rest) = line.strip_prefix("version") {
                if let Some(rest) = rest.trim().strip_prefix('=') {
                    let value = rest.trim().trim_matches('"');
                    if !value.is_empty() {
                        return value.to_string();
                    }
                }
            }
        }
    }
    panic!("no package.version in Cargo.toml");
}

/// Every `src/**/*.rs` file, so the doc checks below cover modules added later without being
/// remembered.
fn sources() -> Vec<PathBuf> {
    let mut out = Vec::new();
    let mut stack = vec![root().join("src")];
    while let Some(dir) = stack.pop() {
        for entry in
            fs::read_dir(&dir).unwrap_or_else(|e| panic!("read_dir {}: {e}", dir.display()))
        {
            let path = entry.expect("dir entry").path();
            if path.is_dir() {
                stack.push(path);
            } else if path.extension().is_some_and(|ext| ext == "rs") {
                out.push(path);
            }
        }
    }
    out.sort();
    assert!(
        out.len() > 10,
        "expected the crate's modules under src/, found {}",
        out.len()
    );
    out
}

// ---- source scanning -------------------------------------------------------------

/// Leading whitespace of `line`, used to tell an `impl` header (column 0) from the methods inside
/// it (indented) without needing to track brace balance across the file.
fn indent_of(line: &str) -> usize {
    line.len() - line.trim_start().len()
}

/// Index just past `word` in `chars` at `i`, if `word` occurs there as a whole word (not as the
/// prefix of a longer identifier — `impl` must not match `implement`).
fn matches_word(chars: &[char], i: usize, word: &str) -> Option<usize> {
    let w: Vec<char> = word.chars().collect();
    if i + w.len() > chars.len() || chars[i..i + w.len()] != w[..] {
        return None;
    }
    if chars
        .get(i + w.len())
        .is_some_and(|c| c.is_alphanumeric() || *c == '_')
    {
        return None;
    }
    Some(i + w.len())
}

fn skip_ws(chars: &[char], mut i: usize) -> usize {
    while i < chars.len() && chars[i].is_whitespace() {
        i += 1;
    }
    i
}

/// Read a `a::b::C` path at index `i`; returns its **last** segment and the index just past it.
/// The last segment is the name the crate root's `pub use` list refers to.
fn read_path(chars: &[char], i: usize) -> Option<(String, usize)> {
    let mut k = i;
    let mut last = String::new();
    loop {
        let start = k;
        while k < chars.len() && (chars[k].is_alphanumeric() || chars[k] == '_') {
            k += 1;
        }
        if k == start {
            break;
        }
        last = chars[start..k].iter().collect();
        if k + 1 < chars.len() && chars[k] == ':' && chars[k + 1] == ':' {
            k += 2;
            continue;
        }
        break;
    }
    if last.is_empty() {
        None
    } else {
        Some((last, k))
    }
}

/// The type an `impl` header implements, as its last path segment: `impl<'a> HullWhite<'a> {` is
/// `HullWhite`, `impl core::fmt::Debug for HullWhite<'_>` is `HullWhite` too (the *type* is what
/// a caller calls the method on, not the trait).  Takes the trimmed header line.
fn impl_target(header: &str) -> Option<String> {
    let chars: Vec<char> = header.chars().collect();
    let after_impl = matches_word(&chars, 0, "impl")?;
    let mut k = skip_ws(&chars, after_impl);
    if chars.get(k) == Some(&'<') {
        // Walk the generic argument list; `<<` / `>>` nest, so count rather than search.
        let mut angle = 0i32;
        while k < chars.len() {
            match chars[k] {
                '<' => angle += 1,
                '>' => angle -= 1,
                _ => {}
            }
            k += 1;
            if angle == 0 {
                break;
            }
        }
    }
    k = skip_ws(&chars, k);
    let (ty, end) = read_path(&chars, k)?;
    // `impl Trait for Type`: the owner is the right-hand side.
    if let Some(after_for) = matches_word(&chars, skip_ws(&chars, end), "for") {
        if let Some((owner, _)) = read_path(&chars, skip_ws(&chars, after_for)) {
            return Some(owner);
        }
    }
    Some(ty)
}

/// The `impl` block the function at `index` is a method of, if any.
///
/// Walks up to the first line indented *strictly less* than the `fn` line and asks whether that line
/// opens an `impl`; blank lines are transparent.  Any other shallower line (`mod`, a `fn`, top-level
/// code) closes the scope, so a free function gets `None`.  Deliberately shallower-line-at-a-time
/// rather than brace-balanced: the crate's methods are `impl` at column 0 with the signatures at
/// four spaces, and a brace counter would have to reason about string literals and prose comments
/// to get the same answer.
fn enclosing_impl(lines: &[&str], index: usize) -> Option<String> {
    let indent = indent_of(lines[index]);
    if indent == 0 {
        return None;
    }
    let mut j = index;
    while j > 0 {
        j -= 1;
        let line = lines[j];
        if line.trim().is_empty() {
            continue;
        }
        if indent_of(line) >= indent {
            continue;
        }
        return impl_target(line.trim_start());
    }
    None
}

/// Pair each line index with the `///` block attached to the item that starts there.
///
/// The pairing has to survive what actually sits between a doc comment and its item: a single-line
/// `#[must_use = "..."]`, a blank line, and a `#[deprecated(since = "...", note = "...")]` spread
/// over four lines.  It is done as one **downward** pass for that reason.  A symmetric upward walk
/// from the `fn` line meets the closing `)]` of a multi-line attribute first and cannot tell it from
/// the end of whatever is above the item, which is how a fully-documented `#[deprecated]`
/// constructor gets reported as undocumented.
///
/// Going down instead: doc lines accumulate into a pending buffer; a line that opens an attribute
/// (`#[`) suspends the buffer until its brackets balance; a blank line leaves the buffer pending; and
/// the next `fn` line takes the buffer.  Any other code line (a struct, a `let`, a closing brace)
/// consumes the buffer without producing an item, because a doc block above it belongs to that
/// item, not to some function further down.
fn doc_blocks(lines: &[&str]) -> HashMap<usize, Vec<String>> {
    let mut out = HashMap::new();
    let mut pending: Vec<String> = Vec::new();
    // > 0 while inside the argument list of an attribute whose opening `#[` is above this line.
    let mut attr_depth: i32 = 0;
    for (i, raw) in lines.iter().enumerate() {
        let line = raw.trim();
        if attr_depth > 0 {
            attr_depth += bracket_delta(line);
            continue;
        }
        if line.starts_with("///") {
            pending.push(line[3..].trim().to_string());
            continue;
        }
        if line.is_empty() {
            continue;
        }
        if line.starts_with("#[") {
            attr_depth += bracket_delta(line);
            continue;
        }
        if fn_name(line).is_some() {
            if !pending.is_empty() {
                out.insert(i, std::mem::take(&mut pending));
            }
            continue;
        }
        pending.clear();
    }
    out
}

/// `(` + `[` + `{` minus `)` + `]` + `}` on one line, with string and character literals removed so
/// a bracket inside a `note = "... (from_yield) ..."` argument does not count as nesting.
fn bracket_delta(line: &str) -> i32 {
    let mut opens = 0i32;
    let mut closes = 0i32;
    let mut in_string = false;
    let mut in_char = false;
    let mut escaped = false;
    for c in line.chars() {
        if escaped {
            escaped = false;
            continue;
        }
        if c == '\\' {
            escaped = true;
            continue;
        }
        if in_string {
            if c == '"' {
                in_string = false;
            }
            continue;
        }
        if in_char {
            if c == '\'' {
                in_char = false;
            }
            continue;
        }
        match c {
            '"' => in_string = true,
            '\'' => in_char = true,
            '(' | '[' | '{' => opens += 1,
            ')' | ']' | '}' => closes += 1,
            _ => {}
        }
    }
    opens - closes
}

/// One function found by scanning the crate's sources, with the facts the checks below ask about.
struct Item {
    path: PathBuf,
    name: String,
    /// The type whose `impl` block this is a method of; `None` for a free function.
    owner: Option<String>,
    /// First line of the attached doc block, empty if there is none.
    summary: String,
    /// The whole doc block.
    doc: Vec<String>,
    /// The signature line starts with `pub ` — necessary, but not sufficient, for being API:
    /// see [`Item::in_public_api`].
    declared_public: bool,
}

impl Item {
    fn label(&self) -> String {
        format!("{}::{}", self.path.display(), self.name)
    }

    /// Whether this function is reachable as `hull_white::…` — the same set `missing_docs` guards.
    ///
    /// Either it lives in a module declared `pub mod` in `src/lib.rs` (so everything `pub` in it is
    /// exported), or its name (free function) or owning type (method) appears in the crate root's
    /// `pub use` list, which is how `HullWhite`'s pricing methods and `get_coupon_times` get out of
    /// otherwise-private modules.
    fn in_public_api(&self, surface: &Surface) -> bool {
        if !self.declared_public {
            return false;
        }
        if self.in_public_module(surface) {
            return true;
        }
        match &self.owner {
            Some(owner) => surface.reexports.contains(owner),
            None => surface.reexports.contains(&self.name),
        }
    }

    fn in_public_module(&self, surface: &Surface) -> bool {
        let relative = self
            .path
            .strip_prefix(root().join("src"))
            .unwrap_or(&self.path)
            .to_string_lossy()
            .replace('\\', "/");
        // `curves.rs`, `curves/mod.rs`, or anything under `curves/`.
        surface.public_modules.iter().any(|m| {
            relative == format!("{m}.rs")
                || relative == format!("{m}/mod.rs")
                || relative.starts_with(&format!("{m}/"))
        })
    }
}

/// The crate's public surface, read out of `src/lib.rs`.
#[derive(Debug)]
struct Surface {
    /// Names the crate root re-exports (`HullWhite`, `get_coupon_times`, `Solution`, ...), after
    /// applying any `as` alias.
    reexports: HashSet<String>,
    /// Modules declared `pub mod` (`curves`, `error`, `test_support`).
    public_modules: Vec<String>,
}

/// Strip `//` comments from a line, ignoring `//` inside string literals (a docs.rs URL in a
/// `#[doc = "..."]` or a `note = "..."` must not truncate the line).
fn strip_comment(line: &str) -> String {
    let mut out = String::new();
    let mut in_string = false;
    let mut escaped = false;
    let chars: Vec<char> = line.chars().collect();
    let mut i = 0;
    while i < chars.len() {
        let c = chars[i];
        if escaped {
            escaped = false;
            out.push(c);
            i += 1;
            continue;
        }
        if c == '\\' {
            escaped = true;
            out.push(c);
            i += 1;
            continue;
        }
        if in_string {
            if c == '"' {
                in_string = false;
            }
            out.push(c);
            i += 1;
            continue;
        }
        if c == '"' {
            in_string = true;
            out.push(c);
            i += 1;
            continue;
        }
        if c == '/' && chars.get(i + 1) == Some(&'/') {
            break;
        }
        out.push(c);
        i += 1;
    }
    out
}

/// Names a single `pub use` statement contributes to the root: the aliases or last path segments
/// inside `pub use a::b::{C, d as e};`, or the one name of `pub use a::b::C;`.
fn use_names(statement: &str) -> Vec<String> {
    let body = statement
        .trim()
        .trim_start_matches("pub use")
        .trim()
        .trim_end_matches(';')
        .trim();
    let mut names = Vec::new();
    if let Some(open) = body.find('{') {
        let close = body.rfind('}').unwrap_or(body.len());
        for piece in body[open + 1..close].split(',') {
            let piece = piece.trim();
            if piece.is_empty() {
                continue;
            }
            let name = match piece.find(" as ") {
                Some(at) => piece[at + " as ".len()..].trim(),
                None => piece,
            };
            if name != "_" && !name.is_empty() {
                names.push(name.to_string());
            }
        }
    } else {
        let name = match body.find(" as ") {
            Some(at) => body[at + " as ".len()..].trim(),
            None => body.rsplit("::").next().unwrap_or(body).trim(),
        };
        if !name.is_empty() && name != "self" {
            names.push(name.to_string());
        }
    }
    names
}

fn surface() -> Surface {
    let lib = read(&root().join("src/lib.rs"));
    let mut reexports = HashSet::new();
    let mut public_modules = Vec::new();
    // A `pub use a::{ ... };` group spans lines, so join the statements first: with comments
    // stripped and newlines folded, a statement runs from "pub use" to the next ';'.
    let folded: String = lib
        .lines()
        .map(|l| format!("{} ", strip_comment(l)))
        .collect();
    let mut from = 0usize;
    while let Some(i) = folded[from..].find("pub use ").map(|i| from + i) {
        let Some(j) = folded[i..].find(';') else {
            break;
        };
        reexports.extend(use_names(&folded[i..i + j + 1]));
        from = i + j + 1;
    }
    for line in lib.lines() {
        let line = strip_comment(line);
        let line = line.trim();
        if let Some(rest) = line.strip_prefix("pub mod ") {
            if let Some((name, _)) = read_path(&rest.chars().collect::<Vec<char>>(), 0) {
                public_modules.push(name);
            }
        }
    }
    let s = Surface {
        reexports,
        public_modules,
    };
    //Fail loudly if the scan found nothing: a silent empty surface would make the public-API checks
    //pass vacuously, which is the failure mode this file exists to prevent.
    assert!(
        !s.reexports.is_empty() && !s.public_modules.is_empty(),
        "could not read the public surface out of src/lib.rs: {s:?}"
    );
    s
}

/// Every function in `src`, with its doc block and owning type resolved.
fn items() -> Vec<Item> {
    let mut out = Vec::new();
    for path in sources() {
        let text = read(&path);
        let lines: Vec<&str> = text.lines().collect();
        let blocks = doc_blocks(&lines);
        for (i, line) in lines.iter().enumerate() {
            let Some(name) = fn_name(line) else { continue };
            let doc = blocks.get(&i).cloned().unwrap_or_default();
            let summary = doc.first().cloned().unwrap_or_default();
            out.push(Item {
                path: path.clone(),
                name,
                owner: enclosing_impl(&lines, i),
                summary,
                doc,
                declared_public: line.trim_start().starts_with("pub "),
            });
        }
    }
    out
}

/// The name in a `fn` signature line, if the line starts one.
fn fn_name(line: &str) -> Option<String> {
    let trimmed = line.trim_start();
    let rest = trimmed
        .strip_prefix("pub ")
        .or_else(|| trimmed.strip_prefix("pub(crate) "))
        .or_else(|| trimmed.strip_prefix("pub(super) "))
        .unwrap_or(trimmed);
    let rest = rest.strip_prefix("const ").unwrap_or(rest);
    let rest = rest.strip_prefix("async ").unwrap_or(rest);
    let rest = rest.strip_prefix("unsafe ").unwrap_or(rest);
    let rest = rest.strip_prefix("fn ")?;
    let name: String = rest
        .chars()
        .take_while(|c| c.is_alphanumeric() || *c == '_')
        .collect();
    if name.is_empty() { None } else { Some(name) }
}

/// Shared walker for the shape checks: `file:name` for every documented function in `src` that
/// trips `is_bad`.
fn offenders(is_bad: impl Fn(&Item) -> bool) -> Vec<String> {
    items()
        .iter()
        // Undocumented functions are `missing_docs`' problem, not these checks'.
        .filter(|item| !item.doc.is_empty())
        .filter(|item| is_bad(item))
        .map(|item| item.label())
        .collect()
}

// ---- release sync ---------------------------------------------------------------

#[test]
fn readme_install_snippet_matches_the_crate_version() {
    let version = crate_version();
    let readme = read(&root().join("README.md"));
    let snippet = format!("hull_white = \"{version}\"");
    assert!(
        readme.contains(&snippet),
        "README.md does not offer the current crate version in its install snippet. \
         Expected `{snippet}`, but package.version is {version}. \
         Bump the README with the version (see \"Doc and version bookkeeping\" in README.md)."
    );
}

#[test]
fn readme_docs_rs_links_name_the_crate_version() {
    let version = crate_version();
    let readme = read(&root().join("README.md"));
    const PREFIX: &str = "docs.rs/hull_white/";
    let mut stale = Vec::new();
    let mut checked = 0usize;
    for (offset, _) in readme.match_indices(PREFIX) {
        let after = &readme[offset + PREFIX.len()..];
        let segment: String = after
            .chars()
            .take_while(|c| *c != '/' && !c.is_whitespace())
            .collect();
        //Only version-shaped segments are links.  The README also *describes* the pattern in prose
        //(`docs.rs/hull_white/<version>/`), and that is not a link to pin.
        if !segment.chars().next().is_some_and(|c| c.is_ascii_digit()) {
            continue;
        }
        checked += 1;
        if segment != version {
            stale.push(segment);
        }
    }
    assert!(
        checked > 0,
        "no versioned docs.rs link found in README.md: the check itself has gone stale"
    );
    assert!(
        stale.is_empty(),
        "README.md links docs.rs at {stale:?} but the crate is {version}. \
         Fix the link with the version bump."
    );
}

#[test]
fn changelog_covers_this_release_and_has_room_for_the_next() {
    let version = crate_version();
    let changelog = read(&root().join("CHANGELOG.md"));
    assert!(
        changelog.contains(&format!("## {version}"))
            || changelog.contains(&format!("## [{version}]")),
        "CHANGELOG.md has no section for the current version {version}. \
         Every release gets one; rename the accumulated ## [Unreleased] block to it."
    );
    assert!(
        changelog.contains("## [Unreleased]") || changelog.contains("## Unreleased"),
        "CHANGELOG.md has no `## [Unreleased]` section for the next release. \
         Add an empty one back after renaming the previous block."
    );
}

// ---- doc-comment shape ----------------------------------------------------------

#[test]
fn no_volality_typos_in_src() {
    //Accepted criterion spelled out as a test so it stays true: the word appeared twice, in
    //`src/model.rs`, on the two functions that share a sentence.
    for path in sources() {
        let text = read(&path);
        assert!(
            !text.contains("volality"),
            "{} still contains \"volality\"; the word is \"volatility\"",
            path.display()
        );
    }
}

#[test]
fn a_now_pricer_is_never_documented_as_pricing_at_a_future_date() {
    //`euro_dollar_future_now` used to open with "Returns price of a Euro Dollar Future at some
    //future time".  A `*_now` function prices today, with the state coming from the calibration
    //(`HullWhite::short_rate_now`); promising a future date on it invites a caller to pass a time
    //that the signature no longer takes.  Checked on the summary line only: a `_now` function may
    //legitimately discuss a future date further down (a forward swap starts in the future even when
    //it is priced now).
    let bad = offenders(|item| {
        let lower = item.summary.to_lowercase();
        let mentions_future = lower.contains("future time")
            || lower.contains("future date")
            || lower.contains("some future");
        item.name.ends_with("_now") && mentions_future
    });
    assert!(
        bad.is_empty(),
        "these *_now functions are summarised as pricing at a future date: {bad:?}"
    );
}

#[test]
fn a_swaption_summary_never_names_the_wrong_side() {
    //`european_receiver_swaption_t` opened with "Returns price of a payer swaption".  On an
    //instrument whose whole identity is which side of the trade you are on, a summary that says the
    //opposite of the function name is worse than no summary: the reader has no way to tell which
    //of the two they are looking at.  Body text naming the other side (a parity note, a
    //cross-reference) is fine; the first line is not.
    let bad = offenders(|item| {
        //A name carrying *both* sides (`swaption_payer_receiver_parity_at_now`, a test of the
        //relation between them) makes the rule vacuous: whatever the summary says, one half of the
        //rule fires.  There is no side for such a function to get wrong.
        if item.name.contains("payer") && item.name.contains("receiver") {
            return false;
        }
        let lower = item.summary.to_lowercase();
        (item.name.contains("receiver") && lower.contains("payer"))
            || (item.name.contains("payer") && lower.contains("receiver"))
    });
    assert!(
        bad.is_empty(),
        "these functions' summary lines name the opposite side of the trade to their own name: {bad:?}"
    );
}

#[test]
fn every_public_function_has_a_summary_that_says_something() {
    //`missing_docs` is denied in `[lints.rust]`, so a *missing* comment on the public surface is a
    //build error.  What it cannot catch is the comment that is present and thin: `/// Price.` passes
    //the lint and tells a caller nothing about which price, of what, valued when.  The proxy for
    //"this summary carries a sentence" is length, and 40 characters is about the shortest thing that
    //can name the instrument and the valuation date.  Scoped to the public API — a private module's
    //helpers are documented for their maintainers, and this is not the check that reviews them.
    let surface = surface();
    let thin: Vec<String> = items()
        .iter()
        .filter(|item| item.in_public_api(&surface))
        .filter(|item| item.summary.chars().count() < 40)
        .map(|item| format!("{} (summary {:?})", item.label(), item.summary))
        .collect();
    assert!(
        thin.is_empty(),
        "these public functions have no summary line, or one too short to say anything: {thin:?}"
    );
}

#[test]
fn the_public_surface_scan_is_not_degenerate() {
    //Guard for the checks above: if `surface()` stopped finding the names that `src/lib.rs`
    //re-exports, `in_public_api` would go false for everything and the thin-summary rule would pass
    //by covering nothing.  These are the exports that matter to this crate's shape today.
    let surface = surface();
    for expected in [
        "HullWhite",
        "get_coupon_times",
        "YieldCurve",
        "from_yield",
        "Solution",
        "SolverSettings",
    ] {
        assert!(
            surface.reexports.contains(expected)
                || surface.public_modules.iter().any(|m| m == expected),
            "`{expected}` is no longer seen in the crate root's public surface: the scan went stale"
        );
    }
    //And the private internals stay out of it, so the rule cannot quietly widen back to "everything".
    let internals = items()
        .iter()
        .filter(|item| item.in_public_api(&surface))
        .filter(|item| {
            item.path.to_string_lossy().contains("/mc.rs")
                || item.path.to_string_lossy().contains("/validation.rs")
        })
        .map(|item| item.label())
        .collect::<Vec<_>>();
    assert!(
        internals.is_empty(),
        "private-module internals ({:?}) are being treated as public API: module visibility changed",
        internals
    );
}
