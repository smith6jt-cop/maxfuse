Save learnings from the current session as a new skill in the skills registry.

## Instructions
1. Review the conversation to identify key learnings:
   - What was the goal?
   - What approaches worked?
   - What failed and why?
   - What exact parameters or configurations succeeded?

2. Create a new skill directory:
   ```
   .skills_registry/plugins/maxfuse/your-skill-name/
   ├── .claude-plugin/
   │   └── plugin.json
   └── skills/
       └── your-skill-name/
           └── SKILL.md
   ```

3. Write `plugin.json` with specific trigger conditions:
   ```json
   {
     "name": "your-skill-name",
     "version": "1.0.0",
     "description": "Specific trigger conditions: (1) condition one, (2) condition two...",
     "author": {"name": "smith6jt"},
     "skills": "./skills"
   }
   ```

4. Write `SKILL.md` with the experiment template format including:
   - Experiment Overview table (date, goal, environment, status)
   - Context section
   - Verified Workflow section
   - **Failed Attempts table** (most valuable section!)
   - Key Insights section

5. Confirm the skill was created successfully
