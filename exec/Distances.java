package  exec;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.UnrecognizedOptionException;

import  cycle.CyclicityTest;
import  tgl.DistancesProcess;

public class Distances {

    /* -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:- */
    public static void main(String[] args) {
        DistancesProcess process = new DistancesProcess();
        Options opts = createOptions();
        String programUsage = "java -jar Distances.jar [opts]";
        try {
            if (processParameters(process, args, opts, programUsage)) {
                process.run();
            }
        } catch (UnrecognizedOptionException uoe) {
            System.out.println(uoe.getMessage());
            HelpFormatter f = new HelpFormatter();
            f.setWidth(120);
            f.printHelp(programUsage, opts);
        } catch (Exception e) {
            System.out.println("Unknown error:");
            e.printStackTrace();
            HelpFormatter f = new HelpFormatter();
            f.setWidth(120);
            f.printHelp(programUsage, opts);
        }
    }


    /* -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:- */
    @SuppressWarnings("static-access")
    private static Options createOptions() {

        String description = null;
        Options opts = new Options();

        /* -------------------------------------------------------------------------------------- */
        /* Cospeciation probability */
        description = "Cospeciation probability.\n" + "Default value = "
                      + DistancesProcess.DEFAULT_COSPECIATION + "\n-\n";
        Option pc = OptionBuilder.withArgName("value").hasArgs(1).withValueSeparator().withDescription(description).create("pc");
        opts.addOption(pc);

        /* -------------------------------------------------------------------------------------- */
        /* Duplication probability */
        description = "Duplication probability.\n" + "Default value = "
                      + DistancesProcess.DEFAULT_DUPLICATION + "\n-\n";
        Option pd = OptionBuilder.withArgName("value").hasArgs(1).withValueSeparator().withDescription(description).create("pd");
        opts.addOption(pd);

        /* -------------------------------------------------------------------------------------- */
        /* Cospeciation probability */
        description = "Host-switch probability.\n" + "Default value = "
                      + DistancesProcess.DEFAULT_HOSTSWITCH + "\n-\n";
        Option ps = OptionBuilder.withArgName("value").hasArgs(1).withValueSeparator().withDescription(description).create("ps");
        opts.addOption(ps);

        /* -------------------------------------------------------------------------------------- */
        /* Cyclicity test */
        description = "Cyclicity test :\n" + CyclicityTest.STOLZER + " - Stolzer et al., 2012\n"
                      + CyclicityTest.DONATI + " - Donati et al., 2014\n" + CyclicityTest.TOFIGH
                      + " - Tofigh et al., 2011\n" + CyclicityTest.TRANSFER_EDGES
                      + " - Only transfer edges\n" + "Default cyclicity test = "
                      + DistancesProcess.DEFAULT_CYCLICITY_TEST + "\n-\n";
        Option test = OptionBuilder.withArgName("value").hasArgs(1).withValueSeparator().withDescription(description).create("ct");
        opts.addOption(test);

        /* -------------------------------------------------------------------------------------- */
        /* Help */
        opts.addOption("h", false, "Print the help message.\n-\n");

        /* -------------------------------------------------------------------------------------- */
        /* Host */
        description = "Host newick file.\n-\n";
        Option host = OptionBuilder.withArgName("file").hasArgs(1).withValueSeparator().withDescription(description).create("host");
        opts.addOption(host);

        /* -------------------------------------------------------------------------------------- */
        /* Host */
        description = "Parasite newick file.\n-\n";
        Option parasite = OptionBuilder.withArgName("file").hasArgs(1).withValueSeparator().withDescription(description).create("parasite");
        opts.addOption(parasite);

        /* -------------------------------------------------------------------------------------- */
        /* Population number of trees */
        description = "Maximum number of trees to generate.\n" + "Default value = "
                      + DistancesProcess.DEFAULT_NUMBER_OF_TREES + "\n-\n";
        Option n = OptionBuilder.withArgName("number").hasArgs(1).withValueSeparator().withDescription(description).create("n");
        opts.addOption(n);

        /* -------------------------------------------------------------------------------------- */
        /* Maximum size factor */
        description = "Maximum size.\n" + "Default value = "
                      + DistancesProcess.DEFAULT_MAXIMUM_SIZE + "\n-\n";
        Option max = OptionBuilder.withArgName("value").hasArgs(1).withValueSeparator().withDescription(description).create("max");
        opts.addOption(max);

        /* -------------------------------------------------------------------------------------- */
        /* Root mapping probability */
        description = "Root mapping probability. Real number in the interval [0.5,1.0].\n"
                      + ("Default value = " + DistancesProcess.DEFAULT_ROOT_MAPPING_PROBABILITY + "\n-\n");
        Option root = OptionBuilder.withArgName("value").hasArgs(1).withValueSeparator().withDescription(description).create("root");
        opts.addOption(root);

        return opts;
    }


    private static boolean processParameters(DistancesProcess process,
                                             String[] args,
                                             Options opts,
                                             String usageString) throws UnrecognizedOptionException,
                                                                ParseException {

        BasicParser bp = new BasicParser();
        CommandLine cl = bp.parse(opts, args);
        if (cl.hasOption("h")) {
            HelpFormatter f = new HelpFormatter();
            f.setWidth(120);
            f.printHelp(usageString, opts);
            return false;
        }

        if (cl.hasOption("pc")) {
            double value = Double.parseDouble(cl.getOptionValue("pc"));
            if (value < 0) {
                throw new UnrecognizedOptionException("Invalid value for option -pc");
            }
            process.setCospeciation(value);
        }

        /* Get duplication probability */
        if (cl.hasOption("pd")) {
            double value = Double.parseDouble(cl.getOptionValue("pd"));
            if (value < 0) {
                throw new UnrecognizedOptionException("Invalid value for option -pd");
            }
            process.setDuplication(value);
        }

        /* Get host-switch probability */
        if (cl.hasOption("ps")) {
            double value = Double.parseDouble(cl.getOptionValue("ps"));
            if (value < 0) {
                throw new UnrecognizedOptionException("Invalid value for option -ps");
            }
            process.setHostSwitch(value);
        }

        /* Get cyclicity test */
        int cyclicityTest = DistancesProcess.DEFAULT_CYCLICITY_TEST;
        if (cl.hasOption("ct")) {
            cyclicityTest = Integer.parseInt(cl.getOptionValue("ct"));
            if (cyclicityTest != CyclicityTest.STOLZER && cyclicityTest != CyclicityTest.TOFIGH
                && cyclicityTest != CyclicityTest.TRANSFER_EDGES
                && cyclicityTest != CyclicityTest.DONATI) {
                throw new UnrecognizedOptionException("Invalid value for option -ct");
            }
            process.setCyclicityTest(cyclicityTest);
        }

        /* Get host file name */
        if (cl.hasOption("host")) {
            process.setHostFile(cl.getOptionValue("host"));
        } else {
            throw new UnrecognizedOptionException("Missing host newick file.");
        }

        /* Get parasite file name */
        if (cl.hasOption("parasite")) {
            process.setParasiteFile(cl.getOptionValue("parasite"));
        } else {
            throw new UnrecognizedOptionException("Missing parasite newick file.");
        }

        /* Get number of trees */
        if (cl.hasOption("n")) {
            int value = Integer.parseInt(cl.getOptionValue("n"));
            if (value <= 0) {
                throw new UnrecognizedOptionException("Invalid value for option -n");
            }
            process.setNumberOfTrees(value);
        }

        /* Get maximum size factor */
        if (cl.hasOption("max")) {
            int value = Integer.parseInt(cl.getOptionValue("max"));
            if (value < 1) {
                throw new UnrecognizedOptionException("Invalid value for option -max");
            }
            process.setMaximumSize(value);
        }

        if (process.getEventProbabilitiesSum() > 1.0) {
            throw new UnrecognizedOptionException("The sum of the probabilities (pc + pd + ps + pl) must be equal to 1.0");
        }

        /* Get root mapping probability */
        if (cl.hasOption("root")) {
            double probability = Double.parseDouble(cl.getOptionValue("root"));
            if (probability < 0.5 || probability > 1.0) {
                throw new UnrecognizedOptionException("Invalid value for option -root");
            }
            process.setRootMappingProbability(probability);
        }

        return true;
    }

}
