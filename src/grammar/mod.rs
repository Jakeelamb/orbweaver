pub mod rule;
pub mod symbol;
pub mod engine;
pub mod digram;
pub mod digram_table;
pub mod builder;

// Export core types for anyone importing the grammar module
pub use crate::grammar::rule::Rule;
pub use crate::grammar::symbol::{Symbol, Direction, SymbolType};
pub use crate::grammar::engine::Grammar;
pub use crate::grammar::builder::GrammarBuilder;
pub use crate::grammar::digram::Digram;
pub use crate::grammar::digram_table::DigramTable;